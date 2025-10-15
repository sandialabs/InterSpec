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
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "InterSpec/PeakFit_imp.hpp"
#include "InterSpec/PeakDists_imp.hpp"
#include "InterSpec/RelActCalcAuto_imp.hpp"

// Undefine isnan and isinf macros
#undef isnan
#undef isinf


using namespace std;

namespace
{
/** To make the code of `PeakFitDiffCostFunction` work with either `PeakDef`
 or `RelActCalcAuto::PeakDefImp<T>`, we'll add in a few interface functions
 */
template<typename T>
struct PeakContinuumImp : public RelActCalcAuto::PeakContinuumImp<T>
{
  std::shared_ptr<const SpecUtils::Measurement> m_external_continuum;
  std::shared_ptr<const SpecUtils::Measurement> externalContinuum() const { return m_external_continuum; }
  void setExternalContinuum( const std::shared_ptr<const SpecUtils::Measurement> &data ){ m_external_continuum = data; }
};

template<typename T>
struct PeakDefImpWithCont : public RelActCalcAuto::PeakDefImp<T>
{
  std::shared_ptr<PeakContinuumImp<T>> m_continuum;
  std::shared_ptr<PeakContinuumImp<T>> getContinuum(){ return m_continuum; }
  void setContinuum( const std::shared_ptr<PeakContinuumImp<T>> &continuum ) { m_continuum = continuum;}
  void setMeanUncert( const T &mean_uncert ){ }
  void setAmplitudeUncert( const T &amp_uncert ){ }
  void setSigmaUncert( const T &sigma_uncert ){ }
  void inheritUserSelectedOptions( const PeakDef &parent, const bool inheritNonFitForValues [[maybe_unused]] )
  {
    //We dont really need this function, since we never use any of this info, but whatever.
    this->m_parent_nuclide = parent.parentNuclide();
    this->m_transition = parent.nuclearTransition();
    this->m_rad_particle_index = parent.decayParticleIndex();
    this->m_xray_element = parent.xrayElement();
    this->m_reaction = parent.reaction();
    this->m_src_energy = (this->m_parent_nuclide || this->m_xray_element || this->m_reaction)
                          ? parent.gammaParticleEnergy() : 0.0f;
    this->m_gamma_type = parent.sourceGammaType();
    //m_rel_eff_index = ...
  }

  static bool lessThanByMean( const PeakDefImpWithCont<T> &lhs, const PeakDefImpWithCont<T> &rhs )
  {
    return (lhs.m_mean < rhs.m_mean);
  }
};//struct PeakDefImpWithCont


/** Replaces contents of input peaks vector with completely new peaks, that crucually have had new continuums allocated. */
void local_unique_copy_continuum( vector<shared_ptr<const PeakDef>> &input_peaks )
{
  map<std::shared_ptr<const PeakContinuum>,vector<shared_ptr<PeakDef>>> contToPeaks;
  for( auto &p : input_peaks )
    contToPeaks[p->continuum()].push_back( make_shared<PeakDef>(*p) );

  for( auto &pp : contToPeaks )
  {
    pp.second[0]->makeUniqueNewContinuum();
    auto newcont = pp.second[0]->continuum();
    for( size_t i = 1; i < pp.second.size(); ++i )
      pp.second[i]->setContinuum( newcont );
  }

  input_peaks.clear();
  for( auto &pp : contToPeaks )
  {
    for( auto p : pp.second )
      input_peaks.push_back( p );
  }
  std::sort( begin(input_peaks), end(input_peaks), &PeakDef::lessThanByMeanShrdPtr );
}//unique_copy_continuum(...)
}//namespace


namespace PeakFitLM
{

/** PeakFitDiffCostFunction is a not exactly apples-to-apples experiment to see how well Ceres does fitting peaks, relative to Minuit2
    This class is somewhat equivalent of the MultiPeakFitChi2Fcn class.
    Or maybe more directly analogous to LinearProblemSubSolveChi2Fcn

 TODO items:
 - [ ] Should add choice of a free-for-all for peak widths, or to have it as a function of energy, like now.
 - [ ] Add in option for a "refine" fit, or an agressive fit
 - [ ] Add in choice to punish for peak being statistically insignificant, like the PeakFitChi2Fcn class
 - [ ] Implement, and check using External continuum, and no-continuum
 - [x] Need to rescale paramters to be around ~1.0.
   - [x] For peak means, should scale from lower to upper roi somehow-ish
   - [x] For FWHM, just scale off initial estimate
 - [x] When using LLS to fit continuum parmaters, shouldnt have these in the fit paramaters.
 - [ ] currently will uses the entire first/last bin or ROI energy range, no matter how-little the ROI extends into those channels - should probably do some rounding or something.
 - [x] Have this class to all its own initial paramater setup, define lower/upper limits, and setup problem; this will make re-using this class a bit easier.
 */
struct PeakFitDiffCostFunction
{
  /** A struct with the info to feed into Ceres. */
  struct ProblemSetup
  {
    /** The starting parameters to use. */
    vector<double> m_parameters;
    /** The indexes of paramters that need to be held constant in the fit.  E.g. if a mean is held constant for a peak whose FWHM/amp is being fit. */
    vector<int> m_constant_parameters;
    /** Lower bounds parameters should be allowed to go to; will not have a value if it shouldnt be restricted. Will be same size as `m_parameters`. */
    vector<std::optional<double>> m_lower_bounds;
    /** Upper bounds parameters should be allowed to go to; will not have a value if it shouldnt be restricted. Will be same size as `m_parameters`. */
    vector<std::optional<double>> m_upper_bounds;
  };//struct ProblemSetup


  PeakFitDiffCostFunction(const std::shared_ptr<const SpecUtils::Measurement> data,
                          const std::vector< std::shared_ptr<const PeakDef> > &starting_peaks,
                          const double roi_lower_energy,
                          const double roi_upper_energy,
                          const PeakContinuum::OffsetType offset_type,
                          const double continuum_ref_energy,
                          const PeakDef::SkewType skew_type,
                          const bool isHPGe,
                          const Wt::WFlags<PeakFitLM::PeakFitLMOptions> options )
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
    m_num_fit_sigmas( ([&starting_peaks]() -> size_t{
      size_t n_fit_sigma = 0;
      for( const auto &p : starting_peaks )
        n_fit_sigma += p->fitFor(PeakDef::Sigma);
      return n_fit_sigma;
    })() ),
    m_use_lls_for_cont( ([&starting_peaks,offset_type]() -> bool {
      if( (offset_type == PeakContinuum::OffsetType::NoOffset) || (offset_type == PeakContinuum::OffsetType::External) )
        return true;
      auto continuum = starting_peaks.empty() ? nullptr : starting_peaks.front()->continuum();
      assert( continuum && (continuum->type() == offset_type) );
      const vector<bool> par_fit_for = continuum ? continuum->fitForParameter() : vector<bool>{};
      for( const bool do_fit : par_fit_for )
      {
        if( !do_fit )
          return false;
      }
      return true;
    })() ),
    m_num_parameters(
      PeakDef::num_skew_parameters(skew_type)
        + (m_use_lls_for_cont ? size_t(0) : PeakContinuum::num_parameters(offset_type))
        + ( options.testFlag(PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent) ? m_num_fit_sigmas : std::min( m_num_fit_sigmas, size_t(2) ) )
        + m_num_peaks
    ),
    m_num_residuals( ((size_t(1) + m_upper_channel) - m_lower_channel)
      + ((!options.testFlag(PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose)
          || options.testFlag(PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig))
          ? (m_num_peaks - 1) : size_t(0))
      + (options.testFlag(PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig) ? size_t(1) : size_t(0))
    ),
    m_offset_type( offset_type ),
    m_ref_energy( continuum_ref_energy ),
    m_skew_type( skew_type ),
    m_external_continuum( ([&starting_peaks,offset_type]() -> shared_ptr<const SpecUtils::Measurement> {
      if( offset_type != PeakContinuum::OffsetType::External )
        return nullptr;
      for( const shared_ptr<const PeakDef> &p : starting_peaks )
      {
        if( (p->continuum()->type() == PeakContinuum::OffsetType::External) && p->continuum()->externalContinuum() )
          return p->continuum()->externalContinuum();
      }
      return nullptr;
    })() ),
    m_max_initial_sigma( ([&starting_peaks](){
      double max_sigma = -1.0;
      for( const shared_ptr<const PeakDef> &p : starting_peaks )
        max_sigma = std::max( max_sigma, (!isinf(p->sigma()) && !isnan(p->sigma())) ? p->sigma() : 1.0 );
      return (max_sigma > 0.01) ? max_sigma : 1.0;
    })() ),
    m_isHPGe( isHPGe ),
    m_options( options ),
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

    const shared_ptr<const vector<float>> &energies = m_data->channel_energies();
    if( !energies || (m_upper_channel > energies->size()) )
      throw runtime_error( "PeakFitDiffCostFunction: m_upper_channel >= energies.size()" );

    if( !m_num_peaks )
      throw runtime_error( "PeakFitDiffCostFunction: no peaks to fit" );
    
    for( size_t i = 1; i < starting_peaks.size(); ++i )
    {
      if( starting_peaks[i]->continuum() != starting_peaks[0]->continuum() )
        throw runtime_error( "PeakFitDiffCostFunction: input peaks must all share a continuum" );
    }

    if( starting_peaks[0]->continuum()->type() != m_offset_type )
      throw runtime_error( "PeakFitDiffCostFunction: input peaks have the same OffsetType as you are fitting for." );

    if( (m_offset_type == PeakContinuum::OffsetType::External)
       && (!m_external_continuum
           || (m_external_continuum->num_gamma_channels() < 7)
           || !m_external_continuum->energy_calibration()
           || !m_external_continuum->energy_calibration()->valid()) )
      throw runtime_error( "PeakFitDiffCostFunction: external continuum wanted, but not provided with starting peaks" );


    auto continuum = m_starting_peaks[0]->continuum();
    // TODO: get rid of m_roi_lower_energy, m_roi_upper_energy, m_offset_type, and m_ref_energy - since they are all duplicated
    assert( (continuum->lowerEnergy() < 5.0)
           || (fabs(continuum->lowerEnergy() - m_roi_lower_energy) < 1.0E-6*std::max(continuum->lowerEnergy(),m_roi_lower_energy)) );
    assert( ((continuum->upperEnergy()+5) > m_data->gamma_energy_max())
           || (fabs(continuum->upperEnergy() - m_roi_upper_energy) < 1.0E-6*std::max(continuum->upperEnergy(),m_roi_upper_energy)) );
    assert( continuum->type() == m_offset_type );
    //assert( continuum->referenceEnergy() == m_ref_energy ); //When dragging edge of ROI to make it smaller, we have to shift the ref enerergy to keep valid
    
    
    
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
    const double num_skew = PeakDef::num_skew_parameters( m_skew_type );
    const size_t num_fit_continuum_pars = m_use_lls_for_cont ? size_t(0) : PeakContinuum::num_parameters( m_offset_type );
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
    if( !m_use_lls_for_cont )
    {
      for( const bool do_fit_par : m_starting_peaks[0]->continuum()->fitForParameter() )
        num_fixed_pars += !do_fit_par;
    }

    // TODO: MultiPeakFitChi2Fcn::dof() adds 1 to the degrees of freedom - not quite sure why that is, at the moment, or if we should do this here.
    return num_channels - 2.0*m_num_peaks + num_fixed_pars - num_sigmas_fit - num_fit_continuum_pars;
  }//double dof() const
  
  size_t number_parameters() const
  {
    assert( m_num_parameters == (PeakDef::num_skew_parameters(m_skew_type)
     + (m_use_lls_for_cont ? size_t(0) : PeakContinuum::num_parameters(m_offset_type))
     + number_sigma_parameters()
     + m_num_peaks) );

    return m_num_parameters;
  }
  
  size_t number_sigma_parameters() const
  {
    if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent) )
      return m_num_fit_sigmas;
    else
      return std::min( m_num_fit_sigmas, size_t(2) );
  }
  
  size_t number_residuals() const
  {
    // If we are punishing for peaks being too close, we will have a `m_num_peaks-1` residuals for
    //  punishing peaks being too close.
    // If we are punishing for insignificant peaks, we will have `m_num_peaks` residuals for this,
    //  but will share the residuals with "being too close" residuals (but then we will need one more).
    assert( m_num_residuals == (
      (((size_t(1) + m_upper_channel) - m_lower_channel)
      + ((!m_options.testFlag(PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose)
          || m_options.testFlag(PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig))
         ? (m_num_peaks - 1) : size_t(0))
      + (m_options.testFlag(PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig) ? size_t(1) : size_t(0))))
    );

    return m_num_residuals;
  }//size_t number_residuals() const


  //
  template<typename T>
  static std::shared_ptr<T> create_continuum( const std::shared_ptr<T> &other_cont [[maybe_unused]]) {
      return std::make_shared<T>();
  }

  template<typename PeakType,typename T>
  vector<PeakType> parametersToPeaks( const T * const params, const T * const uncertainties, T *residuals ) const
  {
#define DEBUG_PAR_TO_PEAKS 0

    const size_t num_fit_continuum_pars = m_use_lls_for_cont ? size_t(0) : PeakContinuum::num_parameters( m_offset_type );
    const size_t end_data_residual = m_upper_channel - m_lower_channel;

    auto starting_continuum = m_starting_peaks[0]->continuum();
    assert( (fabs(starting_continuum->lowerEnergy() - m_roi_lower_energy) < 1.0E-6*std::max(starting_continuum->lowerEnergy(),m_roi_lower_energy))
           || (m_roi_lower_energy < 10.0) );
    assert( (fabs(starting_continuum->upperEnergy() - m_roi_upper_energy) < 1.0E-6*std::max(starting_continuum->upperEnergy(),m_roi_upper_energy))
           || !m_data
           || ((m_roi_upper_energy + 10) > m_data->gamma_energy_max()));
    assert( starting_continuum->type() == m_offset_type );
    //assert( starting_continuum->referenceEnergy() == m_ref_energy ); //may not be true if we are dragging ROI edges to make smaller, so we had to change ref energy to keep it valid


    std::unique_ptr<vector<T>> local_residuals;
    if( !residuals )
    {
      local_residuals.reset( new vector<T>( number_residuals(), T(0.0) ) );
      residuals = local_residuals->data();
    }else
    {
      // Not sure if we need to do this or not
      const size_t nresids = number_residuals();
      for( size_t i = 0; i < nresids; ++i )
        residuals[i] = T(0.0);
    }

    const size_t num_sigmas_fit = number_sigma_parameters();
    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );

    vector<PeakType> peaks, fixed_amp_peaks;
    const double range = m_roi_upper_energy - m_roi_lower_energy;

#if( DEBUG_PAR_TO_PEAKS )
    if constexpr ( !std::is_same_v<T, double> )
    {
      cout << "\nparametersToPeaks: x={";
      for( size_t i = 0; i < number_parameters(); ++i )
        cout << params[i].a << ", ";
      cout << "}" << endl;
    }
#endif

    size_t fit_sigma_num = 0;

    // TODO: if move to fitting the distance between peaks, can remove this next little loop.
    T min_mean = T(std::numeric_limits<double>::infinity());
    T max_mean = T(-std::numeric_limits<double>::infinity());
    for( size_t i = 0; i < m_num_peaks; ++i )
    {
      const size_t mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + i;
      const T frac_roi_range = params[mean_par_index] - 0.5;
      const T mean = m_roi_lower_energy + frac_roi_range*range;
      min_mean = min( min_mean, mean );
      max_mean = max( max_mean, mean );
    }

    for( size_t i = 0; i < m_num_peaks; ++i )
    {
      const size_t mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + i;
      assert( mean_par_index < number_parameters() );

      const bool fit_amp = m_starting_peaks[i]->fitFor(PeakDef::GaussAmplitude);

      // TODO: should move from fitting `frac_roi_range` to instead fitting for the difference between each peak.  Use npeaks unconstrained paramaters, and have thier distances be `exp(par[i])`, then normalize the spacing and add that fraction to each on.  Actually could use npeaks-1 paramaters, and fix the first one to be 0 (so `exp(0) = 1`) to set scale for the other parameters
      const T frac_roi_range = params[mean_par_index] - 0.5;
      const T mean = m_roi_lower_energy + frac_roi_range*range;
      const T amp = fit_amp ? T(1.0) : T(m_starting_peaks[i]->amplitude());

      T sigma, sigma_uncert;

      if( !m_starting_peaks[i]->fitFor(PeakDef::Sigma) || (num_sigmas_fit == 0) )
      {
        sigma = T(m_starting_peaks[i]->sigma());
        sigma_uncert = T(m_starting_peaks[i]->sigmaUncert());
      }else if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent) )
      {
        assert( fit_sigma_num < m_num_fit_sigmas );
        const size_t sigma_index = num_fit_continuum_pars + num_skew + fit_sigma_num;
        assert( sigma_index < number_parameters() );

        sigma = params[sigma_index] * m_max_initial_sigma;
        sigma_uncert = uncertainties ? (uncertainties[sigma_index] * m_max_initial_sigma) : T(0.0);


        if( isnan(sigma) || isinf(sigma) )
        {
          if constexpr ( !std::is_same_v<T, double> )
          {
            cerr << "Recieved a NaN or Inf sigma; params[" << sigma_index<< "] = " << params[sigma_index].a
            << ", sigma=" << sigma.a << endl;
            cerr << "Pars=[";
            for( size_t i = 0; i < m_num_parameters; ++i )
              cerr << params[i].a << ",";
            cerr << "]" << endl;
          }else
          {
            cerr << "Recieved a NaN or Inf sigma; params[" << sigma_index<< "] = " << params[sigma_index]
            << ", sigma=" << sigma << endl;
            cerr << "Pars=[";
            for( size_t i = 0; i < m_num_parameters; ++i )
              cerr << params[i] << ",";
            cerr << "]" << endl;
          }

          throw runtime_error( "Inf or NaN sigma" );
        }//if( isnan(sigma) || isinf(sigma) )
      }else
      {
        // If we fit for only one sigma, we use the paramaters value.
        //  If more than one sigma, then we will start with the paramater value, and adjust
        //  according to the fraction of the ROI the peak is at
        const size_t sigma_index = num_fit_continuum_pars + num_skew;
        sigma = params[sigma_index] * m_max_initial_sigma;
        sigma_uncert = uncertainties ? (uncertainties[sigma_index] * m_max_initial_sigma) : T(0.0);

        if( (i > 0) && (num_sigmas_fit > 1) )
        {
          const size_t first_mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + 0;
          const size_t last_mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + m_num_peaks - 1;
          const T frac_first_roi_range = params[first_mean_par_index] - 0.5;
          const T frac_last_roi_range = params[last_mean_par_index] - 0.5;
          const T first_mean = m_roi_lower_energy + frac_first_roi_range*range;
          const T last_mean = m_roi_lower_energy + frac_last_roi_range*range;

          // I dont know what/why in the world I was doign the below
          T frac_dist;
          //if( (last_mean - first_mean) > 1.0E-2 )
          //  frac_dist = (mean - first_mean) / (last_mean - first_mean);
          //else
          //  frac_dist = T(1.0*i) / (num_sigmas_fit - 1.0);
          assert( mean >= min_mean );
          assert( mean <= max_mean );
          T mean_dists = max_mean - min_mean;
          if( mean_dists > 1.0E-6 )
            frac_dist = (mean - min_mean) / (max_mean - min_mean);
          else
            frac_dist = T(0.5);

          const T fist_sigma = sigma;
          const T last_sigma = sigma * params[sigma_index + 1];
          sigma = fist_sigma + frac_dist*(last_sigma - fist_sigma);

          assert( (sigma_index + 1) < number_parameters() );

          if( uncertainties )
          {
            //TODO: I'm not conviced uncert is correct (note: assuming 100% correlation)
            //TODO: use the covariance between these two parameters to better assign sigma uncert??
            const T first_uncert = sigma_uncert;
            const T last_uncert = sigma_uncert * uncertainties[sigma_index + 1];
            sigma_uncert = first_uncert + frac_dist*(last_uncert - first_uncert);
          }
        }//if( num_sigmas_fit > 1 )

        if( isnan(sigma) || isinf(sigma) )
        {
          if constexpr ( !std::is_same_v<T, double> )
          {
            cerr << "Recieved a NaN _or_ Inf sigma; params[" << sigma_index<< "] = " << params[sigma_index].a
            << ", sigma=" << sigma.a
            << ", num_sigmas_fit=" << num_sigmas_fit
            << endl;
            if( (i > 0) && (num_sigmas_fit > 1) )
            {
              const size_t first_mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + 0;
              const size_t last_mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + m_num_peaks - 1;
              const T frac_first_roi_range = params[first_mean_par_index] - 0.5;
              const T frac_last_roi_range = params[last_mean_par_index] - 0.5;
              const T first_mean = m_roi_lower_energy + frac_first_roi_range*range;
              const T last_mean = m_roi_lower_energy + frac_last_roi_range*range;
              T frac_dist;
              if( (last_mean - first_mean) > 1.0E-2 )
                frac_dist = (mean - first_mean) / (last_mean - first_mean);
              else
                frac_dist = T(1.0*i) / (num_sigmas_fit - 1.0);

              const T fist_sigma = sigma;
              const T last_sigma = sigma * params[sigma_index + 1];
              sigma = fist_sigma + frac_dist*(last_sigma - fist_sigma);

              cerr << "\ti=" << i << ", num_sigmas_fit=" << num_sigmas_fit
              << ", frac_first_roi_range=" << frac_first_roi_range.a
              << ", frac_last_roi_range=" << frac_last_roi_range.a
              << ", first_mean=" << first_mean.a
              << ", last_mean=" << last_mean.a
              << ", frac_dist=" << frac_dist.a
              << ", fist_sigma=" << fist_sigma.a
              << ", last_sigma=" << last_sigma.a
              << ", sigma=" << sigma.a
              << endl;
            }

            cerr << "Pars=[";
            for( size_t i = 0; i < m_num_parameters; ++i )
              cerr << params[i].a << ",";
            cerr << "]" << endl;
          }else
          {
            cerr << "Recieved a NaN _or_ Inf sigma; params[" << sigma_index<< "] = " << params[sigma_index]
            << ", sigma=" << sigma
            << ", num_sigmas_fit=" << num_sigmas_fit
            << endl;
            if( (i > 0) && (num_sigmas_fit > 1) )
            {
              const size_t first_mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + 0;
              const size_t last_mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + m_num_peaks - 1;
              const T frac_first_roi_range = params[first_mean_par_index] - 0.5;
              const T frac_last_roi_range = params[last_mean_par_index] - 0.5;
              const T first_mean = m_roi_lower_energy + frac_first_roi_range*range;
              const T last_mean = m_roi_lower_energy + frac_last_roi_range*range;
              T frac_dist;
              if( (last_mean - first_mean) > 1.0E-2 )
                frac_dist = (mean - first_mean) / (last_mean - first_mean);
              else
                frac_dist = T(1.0*i) / (num_sigmas_fit - 1);
              
              const T fist_sigma = sigma;
              const T last_sigma = sigma * params[sigma_index + 1];
              sigma = fist_sigma + frac_dist*(last_sigma - fist_sigma);

              cerr << "\ti=" << i << ", num_sigmas_fit=" << num_sigmas_fit
              << ", frac_first_roi_range=" << frac_first_roi_range
              << ", frac_last_roi_range=" << frac_last_roi_range
              << ", first_mean=" << first_mean
              << ", last_mean=" << last_mean
              << ", frac_dist=" << frac_dist
              << ", fist_sigma=" << fist_sigma
              << ", last_sigma=" << last_sigma
              << ", sigma=" << sigma
              << endl;
            }//if( (i > 0) && (num_sigmas_fit > 1) )

            cerr << "Pars=[";
            for( size_t i = 0; i < m_num_parameters; ++i )
              cerr << params[i] << ",";
            cerr << "]" << endl;
          }

          throw runtime_error( "Inf or NaN sigma" );
        }//if( isnan(sigma) || isinf(sigma) )

      }//if( num_sigmas_fit == 0 ) / else

      // We will protect against 0 or negative sigmas - even though we should have limited the possible range of
      //  peak widths.  We will nearly-arbitrarily choose a width of 0.2 keV to start checking for reasonableness
      //  wrt spectrum binning width; a nice HPGe might have a 59 keV peak width of 0.25 keV.
      //  We will require peak width to be wider than 0.42 channels; the minimum number of channels per
      //  gaussian sigma that I could find in a commercial detector is 0.5707 channels per Gaussian-sigma.
      //  For a ROI with multiple peaks, we are limiting the Ceres paramaters to be 0.525 channels per gaussian,
      //  but we will let the upper peak vary by 20% either way, and 0.525*0.8=0.42, we will enforce this as our
      //  absolute narrowest allowed peak.
      //  However, if we arent fitting for the FWHM - let it be whatever the user specified.
      if( (sigma < 0.2) && m_starting_peaks[i]->fitFor(PeakDef::Sigma) && (num_sigmas_fit != 0) )
      {
        float peak_mean;
        if constexpr ( std::is_same_v<T, double> )
          peak_mean = static_cast<float>(mean);
        else
          peak_mean = static_cast<float>(mean.a);

        const float left_val = m_data->gamma_channel_lower( m_lower_channel );
        const float upper_val = m_data->gamma_channel_upper( m_upper_channel );
        const double avrg_chnl_width = (upper_val - left_val) / (1.0 + (m_upper_channel - m_lower_channel));
        const double min_reasonable_sigma = (0.42*0.99)*avrg_chnl_width;

        if( (sigma <= min_reasonable_sigma) || (sigma < 0.00025) ) //The 0.00025 is totally arbitrary
        {
          double peak_sigma;
          if constexpr ( std::is_same_v<T, double> )
            peak_sigma = sigma;
          else
            peak_sigma = sigma.a;

          {
            cout << "num_sigmas_fit=" << num_sigmas_fit << ", i=" << i << endl;
            cout << "m_options.testFlag(PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent)=" << m_options.testFlag(PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent) << endl;

            if( (i > 0) && (num_sigmas_fit > 1) )
            {
              const size_t first_mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + 0;
              const size_t last_mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + m_num_peaks - 1;
              const T frac_first_roi_range = params[first_mean_par_index] - 0.5;
              const T frac_last_roi_range = params[last_mean_par_index] - 0.5;
              const T first_mean = m_roi_lower_energy + frac_first_roi_range*range;
              const T last_mean = m_roi_lower_energy + frac_last_roi_range*range;
              T frac_dist;
              if( (last_mean - first_mean) > 1.0E-2 )
                frac_dist = (mean - first_mean) / (last_mean - first_mean);
              else
                frac_dist = T(1.0*i) / (num_sigmas_fit - 1.0);

              const T fist_sigma = sigma;
              const size_t sigma_index = num_fit_continuum_pars + num_skew;
              const T last_sigma = sigma * params[sigma_index + 1];
              const T this_sigma = fist_sigma + frac_dist*(last_sigma - fist_sigma);
              cout << endl;
            }
          }

          throw runtime_error( "peak sigma (" + std::to_string(peak_sigma) + ") smaller than reasonable ("
                              + std::to_string( min_reasonable_sigma )
                              + ", avrg_chnl_width=" + std::to_string(avrg_chnl_width) + ")" );
        }//if( (sigma <= ((0.42*0.99)*avrg_chnl_width)) || (sigma < 0.00025) )
      }//if( sigma is potentually quite small )
#if( DEBUG_PAR_TO_PEAKS )
    if constexpr ( !std::is_same_v<T, double> )
    {
      cout << "parametersToPeaks: Peak " << i << ": mean=" << mean.a << ", sigma=" << sigma.a << endl;
    }
#endif

      if( m_starting_peaks[i]->fitFor(PeakDef::Sigma) )
        fit_sigma_num += 1;

      PeakType peak;
      peak.setMean( mean );
      peak.setSigma( sigma );
      peak.setAmplitude( amp );

      peak.setSkewType( m_skew_type );
      for( size_t skew_index = 0; skew_index < num_skew; ++skew_index )
      {
        const PeakDef::CoefficientType coeff_num = static_cast<PeakDef::CoefficientType>( PeakDef::CoefficientType::SkewPar0 + static_cast<int>(skew_index) );
        const T skew_coef = params[num_fit_continuum_pars + skew_index];
        peak.set_coefficient( skew_coef, coeff_num );
      }

      if( uncertainties )
      {
        T mean_uncert = uncertainties[mean_par_index];
        const T frac_roi_range_uncert = uncertainties[mean_par_index] - 0.5;

        if( !m_starting_peaks[i]->fitFor(PeakDef::Mean) )
          mean_uncert = T(m_starting_peaks[i]->meanUncert());

        if( mean_uncert > 0.0 )
          peak.setMeanUncert( mean_uncert );

        if( sigma_uncert > 0.0 )
          peak.setSigmaUncert( sigma_uncert );
      }//if( uncertainties )

      if( fit_amp )
      {
        peaks.push_back( peak );
      }else
      {
        peak.setAmplitudeUncert( T(m_starting_peaks[i]->amplitudeUncert()) );
        fixed_amp_peaks.push_back( peak );
      }
    }//for( size_t i = 0; i < m_num_peaks; ++i )

    assert( fit_sigma_num == m_num_fit_sigmas );

    // Sort the peaks to make calculating the punishment for peaks being to close easier
    std::sort( begin(peaks), end(peaks), &PeakType::lessThanByMean );


    const size_t nchannel = (m_upper_channel - m_lower_channel) + 1;
    const shared_ptr<const vector<float>> &energies_ptr = m_data->channel_energies();
    assert( energies_ptr && (energies_ptr->size() > (m_lower_channel + nchannel)) );

    assert( !!m_data->gamma_counts() );
    const vector<float> &counts_vec = *m_data->gamma_counts();
    assert( energies_ptr->size() > (m_lower_channel + nchannel) );

    vector<T> peak_counts( nchannel, T(0.0) );
    const float * const energies = energies_ptr->data() + m_lower_channel;
    const float * const channel_counts = counts_vec.data() + m_lower_channel;

    const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_parameters( m_offset_type ) );
    const bool is_step_continuum = PeakContinuum::is_step_continuum( m_offset_type );

    assert( !peaks.empty() || !fixed_amp_peaks.empty() );
    auto continuum = create_continuum( (!peaks.empty()) ? peaks.front().getContinuum() : fixed_amp_peaks.front().getContinuum() );
    continuum->setRange( T(m_roi_lower_energy), T(m_roi_upper_energy) );
    continuum->setType( m_offset_type );

    assert( (m_offset_type != PeakContinuum::OffsetType::External) || !!m_external_continuum );
    assert( (m_offset_type != PeakContinuum::OffsetType::External) || !m_use_lls_for_cont );

    if( m_offset_type == PeakContinuum::OffsetType::External )
      continuum->setExternalContinuum( m_external_continuum );

    if( !m_use_lls_for_cont )
    {
      continuum->setParameters( T(m_ref_energy), params, nullptr );

      for( size_t peak_index = 0; peak_index < peaks.size(); ++peak_index )
        peaks[peak_index].gauss_integral( energies, &(peak_counts[0]), nchannel );

      PeakDists::offset_integral( *continuum, energies, &(peak_counts[0]), nchannel, m_data );
    }else
    {
      vector<T> means, sigmas;
      for( const auto &peak : peaks )
      {
        means.push_back( peak.mean() );
        sigmas.push_back( peak.sigma() );
      }

      const vector<T> skew_pars( params + num_fit_continuum_pars, params + num_fit_continuum_pars + num_skew );

      vector<T> amplitudes, continuum_coeffs, amp_uncerts, cont_uncerts;

      PeakFit::fit_amp_and_offset_imp(energies, channel_counts, nchannel, num_polynomial_terms, is_step_continuum,
                                          T(m_ref_energy), means, sigmas, fixed_amp_peaks, m_skew_type, skew_pars.data(),
                                          amplitudes, continuum_coeffs, amp_uncerts, cont_uncerts,
                                          &(peak_counts[0]) );

      assert( peaks.size() == amplitudes.size() );
      assert( amp_uncerts.empty() || (amp_uncerts.size() == peaks.size()));
      for( size_t peak_index = 0; peak_index < peaks.size(); ++peak_index )
      {
        peaks[peak_index].setAmplitude( amplitudes[peak_index] );
        if( !amp_uncerts.empty() )
          peaks[peak_index].setAmplitudeUncert( amp_uncerts[peak_index] );
      }

      if( (m_offset_type != PeakContinuum::OffsetType::NoOffset)
         && (m_offset_type != PeakContinuum::OffsetType::External) )
      {
        continuum->setParameters( T(m_ref_energy), continuum_coeffs.data(), cont_uncerts.data() );
      }
    }//if( m_use_lls_for_cont )


    // Put fixed amplitude peaks into `peaks`.
    if( !fixed_amp_peaks.empty() )
    {
      for( auto &peak : fixed_amp_peaks )
        peaks.push_back( peak );
      std::sort( begin(peaks), end(peaks), &PeakType::lessThanByMean );
    }


    for( auto &peak : peaks )
      peak.setContinuum( continuum );


    T chi2( 0.0 );
    for( size_t residual_num = 0; residual_num < nchannel; ++residual_num )
    {
      const double ndata = channel_counts[residual_num];

      assert( residual_num < number_residuals() );

      if( ndata >= PEAK_FIT_MIN_CHANNEL_UNCERT )
        residuals[residual_num] = (ndata - peak_counts[residual_num]) / sqrt(ndata);
      else
        residuals[residual_num] = peak_counts[residual_num]; // TODO: This is ad-hoc, and just following PeakFitChi2Fcn implementation for the moment

      chi2 += residuals[residual_num] * residuals[residual_num];
    }//for( size_t channel = m_lower_channel; channel <= m_upper_channel; ++channel )

    const T chi2_dof = chi2 / dof();
    for( size_t i = 0; i < m_num_peaks; ++i )
      peaks[i].set_coefficient( chi2_dof, PeakDef::Chi2DOF );

    // Punish if the peaks are too close.
    for( size_t i = 1; i < peaks.size(); ++i )
    {
      T sigma;
      if constexpr ( !std::is_same_v<T, double> )
      {
        const T sigma_i = peaks[i].sigma();
        const T sigma_j = peaks[i-1].sigma();
        sigma = 0.5*(sigma_j + sigma_i);
      }else
      {
        const double sigma_i = peaks[i].gausPeak() ? peaks[i].sigma() : 0.25*peaks[i].roiWidth();
        const double sigma_j = peaks[i-1].gausPeak() ? peaks[i-1].sigma() : 0.25*peaks[i-1].roiWidth();
        sigma = 0.5*(sigma_j + sigma_i);
      }


      if( !m_options.testFlag(PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose) )
      {
        assert( i >= 1 );
        const T dist = abs(peaks[i-1].mean() - peaks[i].mean());
        T reldist = dist / sigma;
        if( (reldist < 0.01) || isinf(reldist) || isnan(reldist) )
        {
          if constexpr ( !std::is_same_v<T, double> )
          {
            reldist.a = 0.01;
            reldist.v = peaks[i-1].mean().v - peaks[i].mean().v;
          }else
          {
            reldist = 0.01;
          }
        }//if( (reldist < 0.01) || isinf(reldist) || isnan(reldist) )

        // I'm not really sure what should be used to weight the punishement of peaks being too close...
        const double punishment_factor = 0.25 * (m_upper_channel - m_lower_channel) / m_num_peaks;
        //const double punishment_factor = 5;

        assert( (end_data_residual + i) <= number_residuals() );

        if( reldist < 1.25 )
        {
          residuals[end_data_residual + i - 1] = (punishment_factor / reldist);
        }
        //else if( reldist < 1.5 )
        //{
        //  // At reldist==1, this will be equal to the above case,
        //  //  and by reldist==1.5: exp( -pow(5*0.5,3.0) )== 1.64E-7
        //  residuals[last_data_residual + i] = punishment_factor * exp( -pow(5.0*(reldist - 1.0),3.0) );
        //}
        else
        {
          residuals[end_data_residual + i - 1] = T(0.0);
        }
      }//


      if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig) )
      {
        // As of 20250415 - this section of code is totally untested.
        //From LinearProblemSubSolveChi2Fcn:
        //If the peak area is statistically insignificant on the interval
        //  -1.75sigma to 1.75, then punish!

        // Note: for the `PeakFitLMOptions::PunishForPeakBeingStatInsig` option, we ad a residual
        //       over
        assert( (end_data_residual + peaks.size()) == number_residuals() );
        residuals[end_data_residual + peaks.size() - 1] = T(0.0); //We didnt initialize this residual yet

        for( size_t i = 0; i < peaks.size(); ++i )
        {
          const PeakType &peak = peaks[i];

          if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose) )
            residuals[end_data_residual + i] = T(0.0);

          double mean, sigma, amp;
          if constexpr ( std::is_same_v<T, double> )
          {
            mean = peak.mean();
            sigma = peak.sigma();
            amp = peak.amplitude();
          }else
          {
            mean = peak.mean().a;
            sigma = peak.sigma().a;
            amp = peak.amplitude().a;
          }

          const float lower_energy = static_cast<float>( mean - 1.75*sigma );
          const float upper_energy = static_cast<float>( mean + 1.75*sigma );

          const size_t lower_channel = std::max(m_lower_channel, m_data->find_gamma_channel(lower_energy) );
          const size_t upper_channel = std::min(m_upper_channel, m_data->find_gamma_channel(upper_energy) );

          double dataarea = 0.0;
          for( size_t bin = lower_channel; (bin < counts_vec.size()) && (bin < upper_channel); ++bin )
            dataarea += counts_vec[bin];

          // TODO: all this punishment hasnt been tested, or thought out
          const double punishment_chi2 = (2.0*(1.0 + upper_channel - lower_channel)); //

          if( amp < std::max(2.0*sqrt(dataarea), 1.0) )
          {
            if( amp <= 1.0 )
            {
              if constexpr ( std::is_same_v<T, double> )
              {
                residuals[end_data_residual + i] += 100.0 * punishment_chi2;
              }else
              {
                T punish( 1.0 );
                punish.v = peak.amplitude().v / peak.amplitude().a; // I think this gets the matrix part to ~unity amplitude
                punish *= 100.0 * punishment_chi2;
                residuals[end_data_residual + i] += punish;
              }
            }else
            {
              residuals[end_data_residual + i] += (T(-1.0) + 2.0*sqrt(dataarea)/peak.amplitude()) * punishment_chi2;
            }
          }//if( peaks[i].amplitude() < 2.0*sqrt(dataarea) )
        }//for( size_t i = 0; i < peaks.size(); ++i )
      }//if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig) )
    }//for( size_t i = 1; i < peaks.size(); ++i )


    // We could skip pretty much the rest of the function, if we are using `PeakDefImpWithCont`,
    //  but we'll just leave it in to make sure equal treatment of things.
    //
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
      const T new_mean = peaks[new_index].mean();

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
          if( abs(new_mean - orig_mean) < abs(new_mean - prev_mean) )
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

#if( !defined(NDEBUG) && PERFORM_DEVELOPER_CHECKS && defined(CERES_PUBLIC_JET_H_) )
    const size_t nresids = number_residuals();
    for( size_t i = 0; i < nresids; ++i )
    {
      if constexpr ( std::is_same_v<T, double> )
      {
        assert( !isnan(residuals[i]) && !isinf(residuals[i]) );
      }else
      {
        check_jet_for_NaN( residuals[i] );
      }
    }//for( size_t i = 0; i < nresids; ++i )
#endif

    return peaks;
  }


  template<typename T>
  bool operator()( T const *const *parameters, T *residuals ) const
  {
    m_ncalls += 1;

    try
    {
      T const * const pars = parameters[0];
      const T * const uncertainties = nullptr;
      if constexpr ( std::is_same_v<T, double> )
      {
        parametersToPeaks<PeakDef,double>( pars, uncertainties, residuals );
      }else
      {
        parametersToPeaks<PeakDefImpWithCont<T>,T>( pars, uncertainties, residuals );
      }
    }catch( std::exception &e )
    {
      cerr << "PeakFitDiffCostFunction: caught exception: '" << e.what() << "'" << endl;
      return false;
    }

    return true;
  }//operator()


  static PeakDef::SkewType skew_type_from_prev_peaks( const vector<shared_ptr<const PeakDef> > &inpeaks )
  {
    PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;

    for( const auto &p : inpeaks )
    {
      // PeakDef::SkewType is defined roughly in order of how we should prefer them.
      //  However, if one of the peaks in the ROI required a higher-level of skew function
      //  to describe it, we will use that skew for the entire ROI.
      // Here we will use the skew value of the first peak we encounter, with the higher-level
      //  of skew function, as the starting value, and wether we should fit for the parameter
      //  at all.
      if( p->skewType() > skew_type )
        skew_type = p->skewType();
    }//for( const auto &p : inpeaks )

    return skew_type;
  }//void skew_type_from_prev_peaks(...)


  void setup_skew_parameters( double *pars,
                             vector<int> &constant_parameters,
                             vector<std::optional<double>> &lower_bounds,
                             vector<std::optional<double>> &upper_bounds,
                          const vector<shared_ptr<const PeakDef> > &inpeaks ) const
  {
    const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skew_type );

    if( num_skew_pars == 0 )
      return;

    vector<bool> fit_parameter( num_skew_pars, true );
    vector<double> starting_value( num_skew_pars, 0.0 );
    vector<double> lower_values( num_skew_pars, 0.0 ), upper_values( num_skew_pars, 0.0 );

    for( size_t i = 0; i < num_skew_pars; ++i )
    {
      const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i);

      double lower, upper, start, dx;
      const bool use = PeakDef::skew_parameter_range( m_skew_type, ct, lower, upper, start, dx );
      assert( use );
      if( !use )
        throw logic_error( "Inconsistent skew par val" );

      starting_value[i] = start;
      lower_values[i] = lower;
      upper_values[i] = upper;
    }//for( size_t i = 0; i < num_skew_pars; ++i )


    for( const auto &p : inpeaks )
    {
      if( p->skewType() != m_skew_type )
        continue;

      for( size_t i = 0; i < num_skew_pars; ++i )
      {
        const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i);
        fit_parameter[i] = p->fitFor( ct );
        double val = p->coefficient( ct );

        if( IsInf(val) || IsNan(val) || (val < lower_values[i]) || (val > upper_values[i]) )
          val = starting_value[i];

        // When we are dragging a ROI edge, sometimes the power-law can get stuck up around 100,
        //  so lets try to avoid this
        switch( m_skew_type )
        {
          case PeakDef::NoSkew:   case PeakDef::Bortel:
          case PeakDef::GaussExp: case PeakDef::ExpGaussExp:
            break;

          case PeakDef::CrystalBall:
          case PeakDef::DoubleSidedCrystalBall:
          {
            switch( ct )
            {
              case PeakDef::Mean:           case PeakDef::Sigma:
              case PeakDef::GaussAmplitude: case PeakDef::NumCoefficientTypes:
              case PeakDef::Chi2DOF:

              case PeakDef::SkewPar0:
              case PeakDef::SkewPar2:
                if( (val > 3.0) && fit_parameter[i] )
                  val = starting_value[i];
                break;

              case PeakDef::SkewPar1:
              case PeakDef::SkewPar3:
                if( (val > 6.0) && fit_parameter[i] )
                  val = starting_value[i];
                break;
            }//switch( ct )
            break;
          }//A Crystal Ball distribution
        }//switch( m_skew_type )

        starting_value[i] = val;
      }//for( size_t i = 0; i < num_skew_pars; ++i )
    }//for( const auto &p : inpeaks )

    assert( num_skew_pars == fit_parameter.size() );
    assert( num_skew_pars == starting_value.size() );

    const size_t num_fit_continuum_pars = m_use_lls_for_cont ? size_t(0) : PeakContinuum::num_parameters( m_offset_type );

    for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
    {
      const size_t par_index = num_fit_continuum_pars + skew_index;
      pars[par_index] = starting_value[skew_index];

      if( fit_parameter[skew_index] )
      {
        lower_bounds[par_index] = std::max( lower_values[skew_index], 0.5*starting_value[skew_index] );
        upper_bounds[par_index] = std::min( upper_values[skew_index], 1.5*fabs(starting_value[skew_index]) );
      }else
      {
        constant_parameters.push_back( static_cast<int>(par_index) );
      }
    }//for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
  }//setup_skew_parameters(...)




  ProblemSetup get_problem_setup() const
  {
    const size_t num_fit_pars = number_parameters();
    const size_t num_sigmas_fit = number_sigma_parameters();
    const size_t num_fit_continuum_pars = m_use_lls_for_cont ? size_t(0) : PeakContinuum::num_parameters( m_offset_type );
    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );

    vector<int> constant_parameters;
    vector<double> parameters( num_fit_pars, 0.0 );
    vector<std::optional<double>> lower_bounds( num_fit_pars );
    vector<std::optional<double>> upper_bounds( num_fit_pars );

    if( m_use_lls_for_cont )
    {
      assert( num_fit_continuum_pars == 0 );
    }else
    {
      const shared_ptr<const PeakContinuum> initial_continuum = m_starting_peaks.front()->continuum();

      assert( initial_continuum );
      const vector<double> &cont_pars = initial_continuum->parameters();
      const vector<bool> par_fit_for = initial_continuum->fitForParameter();
      assert( cont_pars.size() == num_fit_continuum_pars );
      assert( par_fit_for.size() == num_fit_continuum_pars );

      for( size_t i = 0; i < num_fit_continuum_pars; ++i )
      {
        parameters[i] = cont_pars[i];
        if( !par_fit_for[i] )
          constant_parameters.push_back( static_cast<int>(i) );
      }
    }//if( cost_function->m_use_lls_for_cont )


    setup_skew_parameters( &(parameters[0]), constant_parameters, lower_bounds, upper_bounds, m_starting_peaks );

    const double range = (m_roi_upper_energy - m_roi_lower_energy);
    const size_t midbin = m_data->find_gamma_channel( 0.5*(m_roi_lower_energy + m_roi_upper_energy) );
    const float binwidth = m_data->gamma_channel_width( midbin );

    const double avrg_bin_width = ([this](){
      const size_t lower_bin = m_data->find_gamma_channel(m_roi_lower_energy);
      const size_t upper_bin = m_data->find_gamma_channel(m_roi_upper_energy);
      const double left_energy = m_data->gamma_channel_lower(lower_bin);
      const double right_energy = m_data->gamma_channel_upper(upper_bin);
      return (right_energy - left_energy) / (1 + (upper_bin - lower_bin));
    })();

    double minsigma = binwidth, min_input_sigma = binwidth;
    double maxsigma = 0.5*range, max_input_sigma = binwidth;


    for( size_t i = 0; i < m_starting_peaks.size(); ++i )
    {
      const std::shared_ptr<const PeakDef> &peak = m_starting_peaks[i];
      const double mean  = peak->mean();
      const double sigma = peak->sigma();

      const size_t mean_par_index = num_fit_continuum_pars + num_skew + num_sigmas_fit + i;

      // TODO:
      // TODO: for more than one peak, move to using parameters to, move to using npeaks unconstrained paramaters, and have thier distances be `exp(par[i])` from the previous peak, then normalize the spacing and add that fraction to each on.  Actually could use npeaks-1 paramaters, and fix the first one to be 0 (so `exp(0) = 1`) to set scale for the other parameters
      const double rel_mean = 0.5 + (mean - m_roi_lower_energy) / range;
      parameters[mean_par_index] = rel_mean;

      if( !peak->fitFor(PeakDef::Mean) )
      {
        constant_parameters.push_back( static_cast<int>(mean_par_index) );
      }else
      {
        // This needs some much better limits

        //        if( coFitPeaks[i] == candidatepeak )
        //        {
        //          lower_bounds[mean_par_index] = rel_mean - sigma/range;
        //          upper_bounds[mean_par_index] = rel_mean + sigma/range;
        //        }else
        //        {
        //          lower_bounds[mean_par_index] = rel_mean - 0.5*sigma/range;
        //          upper_bounds[mean_par_index] = rel_mean + 0.5*sigma/range;
        //        }

        // If we want to limit the peaks to be within ROI, we would limit its range to [0.5, 1.5],
        //  however, we might be fitting a peak outside the ROI, so a little below, if the
        //  peak is already outside the ROI, we'll let it roam around a little there:
        lower_bounds[mean_par_index] = 0.5;
        upper_bounds[mean_par_index] = 1.5;

        if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::MediumRefinementOnly) )
        {
          const double plus_50p_sigma = mean + 0.5*sigma;
          const double minus_50p_sigma = mean - 0.5*sigma;
          lower_bounds[mean_par_index] = 0.5 + (minus_50p_sigma - m_roi_lower_energy) / range;
          upper_bounds[mean_par_index] = 0.5 + (plus_50p_sigma - m_roi_lower_energy) / range;
          assert( rel_mean >= *lower_bounds[mean_par_index] );
          assert( rel_mean <= *upper_bounds[mean_par_index] );
        }

        if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::SmallRefinementOnly) )
        {
          const double plus_15p_sigma = mean + 0.15*sigma;
          const double minus_15p_sigma = mean - 0.15*sigma;
          lower_bounds[mean_par_index] = 0.5 + (minus_15p_sigma - m_roi_lower_energy) / range;
          upper_bounds[mean_par_index] = 0.5 + (plus_15p_sigma - m_roi_lower_energy) / range;
          assert( rel_mean >= *lower_bounds[mean_par_index] );
          assert( rel_mean <= *upper_bounds[mean_par_index] );
        }

        // If the peak is already outside the ROI - we'll exand the limits, and give it a little extra to roam
        if( rel_mean <= *lower_bounds[mean_par_index] )
          lower_bounds[mean_par_index] = rel_mean - std::max( 0.1, fabs(0.25*rel_mean) ); // A little careful to not g onegative
        if( rel_mean >= *upper_bounds[mean_par_index] )
          upper_bounds[mean_par_index] = 1.2*rel_mean;

        assert( rel_mean >= *lower_bounds[mean_par_index] );
        assert( rel_mean <= *upper_bounds[mean_par_index] );
      }//if( fixed mean ) / else

      min_input_sigma = std::min( min_input_sigma, sigma );
      max_input_sigma = std::max( max_input_sigma, sigma );


      if( !peak->fitFor(PeakDef::Sigma) )
      {
        minsigma = std::min( minsigma, sigma );
        maxsigma = std::max( maxsigma, sigma );
      }else
      {
        float lowersigma, uppersigma;
        expected_peak_width_limits( mean, m_isHPGe, m_data, lowersigma, uppersigma );
        if( i == 0 )
          minsigma = lowersigma;
        if( i == (m_starting_peaks.size()-1) )
          maxsigma = uppersigma;
      }
    }//for( size_t i = 0; i < coFitPeaks.size(); ++i )

    if( num_sigmas_fit == 0 )
    {
      // We have no paramaters dedicated to fitting sigmas
    }else if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent) )
    {
      // We have one paramater for every sigma we are fitting
      size_t fit_sigma_index = 0;
      for( size_t i = 0; i < m_starting_peaks.size(); ++i )
      {
        const std::shared_ptr<const PeakDef> &peak = m_starting_peaks[i];
        if( peak->fitFor(PeakDef::Sigma) )
        {
          const double sigma = peak->sigma();
          const size_t sigma_index = num_fit_continuum_pars + num_skew + fit_sigma_index;

          parameters[sigma_index] = sigma / m_max_initial_sigma;
          lower_bounds[sigma_index] = (0.5*std::min(minsigma,min_input_sigma)) / m_max_initial_sigma;
          upper_bounds[sigma_index] = (1.5*std::max(maxsigma,max_input_sigma)) / m_max_initial_sigma;

          if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::MediumRefinementOnly) )
          {
            lower_bounds[sigma_index] = 0.5 * sigma / m_max_initial_sigma;
            upper_bounds[sigma_index] = 1.5 * sigma / m_max_initial_sigma;
          }

          if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::SmallRefinementOnly) )
          {
            lower_bounds[sigma_index] = 0.85 * sigma / m_max_initial_sigma;
            upper_bounds[sigma_index] = 1.15 * sigma / m_max_initial_sigma;
          }

          assert( parameters[sigma_index] >= *lower_bounds[sigma_index] );
          assert( parameters[sigma_index] <= *upper_bounds[sigma_index] );

          fit_sigma_index += 1;
        }//if( peak->fitFor(PeakDef::Sigma) )
      }//for( size_t i = 0; i < m_starting_peaks.size(); ++i )

      assert( fit_sigma_index == num_sigmas_fit );
    }else
    {
      // We have one paramater if fitting a single sigma, or two paramaters if fitting multiple.
      const size_t fit_sigma_index = num_fit_continuum_pars + num_skew;
      parameters[fit_sigma_index] = 1.0;

      lower_bounds[fit_sigma_index] = (0.5*std::min(minsigma,0.75*min_input_sigma)) / m_max_initial_sigma;
      upper_bounds[fit_sigma_index] = (1.5*std::max(maxsigma,1.25*max_input_sigma)) / m_max_initial_sigma;

      assert( parameters[fit_sigma_index] >= *lower_bounds[fit_sigma_index] );
      assert( parameters[fit_sigma_index] <= *upper_bounds[fit_sigma_index] );

      if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::MediumRefinementOnly) )
      {
        // These bounds are only approximately the claimed 50% - since we're lumping all peaks in ROI together
        //  but its probably not worth the hassle of correctness to go into much more detail
        lower_bounds[fit_sigma_index] = (0.5*min_input_sigma) / m_max_initial_sigma;
        upper_bounds[fit_sigma_index] = (1.5*max_input_sigma) / m_max_initial_sigma;

        assert( parameters[fit_sigma_index] >= *lower_bounds[fit_sigma_index] );
        assert( parameters[fit_sigma_index] <= *upper_bounds[fit_sigma_index] );
      }

      if( m_options.testFlag(PeakFitLM::PeakFitLMOptions::SmallRefinementOnly) )
      {
        // These bounds are only approximately the claimed 15% - since we're lumping all peaks in ROI together
        //  but its probably not worth the hassle of correctness to go into much more detail
        lower_bounds[fit_sigma_index] = (0.85*min_input_sigma) / m_max_initial_sigma;
        upper_bounds[fit_sigma_index] = (1.15*max_input_sigma) / m_max_initial_sigma;

        assert( parameters[fit_sigma_index] >= *lower_bounds[fit_sigma_index] );
        assert( parameters[fit_sigma_index] <= *upper_bounds[fit_sigma_index] );
      }

      // The least number of channels per FWHM that I could find in commercial detectors is a FWHM of 1.344 channels,
      //  so we will make this an absolute lower bound on FWHM (we are actually fitting gaussian sigma, so this
      //  becomes 0.5707 channels per sigma - so we'll round down to 0.525, which since peaks can go down to 0.8
      //  times this, brings down the absolute narrowest peak to a sigma of 0.42 channels - narrow than reasonable)
      //
      //  If any peak has a fixed width that is less than `0.525*avrg_bin_width`, then we wont apply this constraint
      bool any_peak_fixed_narrow = false, all_widths_fixed = true;
      for( size_t i = 0; i < m_starting_peaks.size(); ++i )
      {
        const std::shared_ptr<const PeakDef> &peak = m_starting_peaks[i];
        const bool fit_width = peak->fitFor(PeakDef::CoefficientType::Sigma);
        any_peak_fixed_narrow |= (!fit_width && (peak->sigma() < 0.525*avrg_bin_width));
        all_widths_fixed = (all_widths_fixed && !fit_width);
      }

      if( !any_peak_fixed_narrow )
      {
        const double rel_bin_width = avrg_bin_width / m_max_initial_sigma;
        //cout << "Before: par[fit_sigma_index]={" << lower_bounds[fit_sigma_index].value()
        //<< " < " << parameters[fit_sigma_index]
        //<< " < " << upper_bounds[fit_sigma_index].value()
        //<< "}, with 0.525*rel_bin_width=" << 0.525*rel_bin_width << endl;

        lower_bounds[fit_sigma_index] = std::max( 0.525*rel_bin_width, lower_bounds[fit_sigma_index].value() );
        upper_bounds[fit_sigma_index] = std::max( 0.525*rel_bin_width, upper_bounds[fit_sigma_index].value() );
        parameters[fit_sigma_index]   = std::max( 0.525*rel_bin_width, parameters[fit_sigma_index] );
      }//if( !any_peak_fixed_narrow )

      assert( parameters[fit_sigma_index] >= *lower_bounds[fit_sigma_index] );
      assert( parameters[fit_sigma_index] <= *upper_bounds[fit_sigma_index] );
      assert( all_widths_fixed || (lower_bounds[fit_sigma_index].has_value() && (lower_bounds[fit_sigma_index].value() >= 0.525*(avrg_bin_width / m_max_initial_sigma)) ) );

      if( num_sigmas_fit > 1 )
      {
        const size_t upper_sigma_index = fit_sigma_index + 1;

        //Allow FWHM to vary by +-20% over the ROI.
        // The `fit_sigma_index` parameters is the width of the first peak in ROI,
        //  and `upper_sigma_index` is the multiplier of that first ROI, to give
        //  the width of the last peak in the ROI.  Peaks between the first and last peak are scaled
        //  according to thier distance between.
        parameters[upper_sigma_index] = 1.0;
        lower_bounds[upper_sigma_index] = 0.8;
        upper_bounds[upper_sigma_index] = 1.2;

        // TODO: if PeakFitLM::PeakFitLMOptions::MediumRefinementOnly or PeakFitLM::PeakFitLMOptions::SmallRefinementOnly
        //       we should limit things in a bit more sensical way

        assert( parameters[upper_sigma_index] >= *lower_bounds[upper_sigma_index] );
        assert( parameters[upper_sigma_index] <= *upper_bounds[upper_sigma_index] );
      }
    }//if( !num_sigmas_fit ) / else

    ProblemSetup prob_setup;
    prob_setup.m_parameters = parameters;
    prob_setup.m_constant_parameters = constant_parameters;
    prob_setup.m_lower_bounds = lower_bounds;
    prob_setup.m_upper_bounds = upper_bounds;

    return prob_setup;
  }//const ProblemSetup get_problem_setup() const

public:
  const std::shared_ptr<const SpecUtils::Measurement> m_data;
  const size_t m_lower_channel;
  const size_t m_upper_channel;
  const double m_roi_lower_energy;
  const double m_roi_upper_energy;
  
  const size_t m_num_peaks;
  const std::vector< std::shared_ptr<const PeakDef> > m_starting_peaks;
  const size_t m_num_fit_sigmas;
  const bool m_use_lls_for_cont;
  const size_t m_num_parameters;
  const size_t m_num_residuals;
  const PeakContinuum::OffsetType m_offset_type;
  const double m_ref_energy;
  const PeakDef::SkewType m_skew_type;
  const std::shared_ptr<const SpecUtils::Measurement> m_external_continuum;
  const double m_max_initial_sigma;
  const bool m_isHPGe;
  const Wt::WFlags<PeakFitLMOptions> m_options;

  mutable std::atomic<unsigned int> m_ncalls;
};//struct PeakFitDiffCostFunction



/** All peaks passed in must share a PeakContinuum.
 */
vector<shared_ptr<const PeakDef>> fit_peaks_in_roi_LM( const vector<shared_ptr<const PeakDef>> coFitPeaks,
                                                      const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                                      const bool isHPGe,
                                                      const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options )
{
  /** For this first go, we will have Ceres fit for things.

   This is only tested to the "seems to work, for simple cases" level.

   In the future:
     - Need to deal with errors during solving, and during Covariance computation
     - Deal with setting number of threads reasonably.
     - Try out using a ceres::LossFunction
     - Optimize/pcik the Ceres-based paramaters to get reasonable fits.
   */

  if( coFitPeaks.empty() )
    throw runtime_error( "fit_peaks_in_roi_LM: empty input peaks." );

  // Lets make sure all coFitPeaks share a continuum.
  for( size_t i = 1; i < coFitPeaks.size(); ++i )
  {
    if( coFitPeaks[i]->continuum() != coFitPeaks[0]->continuum() )
      throw runtime_error( "fit_peak_for_user_click_LM: input peaks all must share a single continuum" );
  }//for( size_t i = 1; i < coFitPeaks.size(); ++i )

  try
  {
    // TODO: check that repeaded calls to this function wont cause adding/removing channels due to find_gamma_channel(...) rounding type things - e.g., do we need to subtract off a tiny bit from roiUpperEnergy

    const shared_ptr<const PeakContinuum> input_continuum = coFitPeaks[0]->continuum();

    const float roiLowerEnergy = static_cast<float>( input_continuum->lowerEnergy() );
    const float roiUpperEnergy = static_cast<float>( input_continuum->upperEnergy() );
    assert( roiLowerEnergy < roiUpperEnergy );
    if( roiLowerEnergy >= roiUpperEnergy )
      throw runtime_error( "Invalid energy range (" + std::to_string(roiLowerEnergy) + ", " + std::to_string(roiUpperEnergy) + ")" );


    const size_t lower_channel = dataH->find_gamma_channel(roiLowerEnergy);
    const size_t upper_channel = dataH->find_gamma_channel(roiUpperEnergy);
    const double roi_lower = dataH->gamma_channel_lower( lower_channel );
    const double roi_upper = dataH->gamma_channel_upper( upper_channel );
    
    assert( (lower_channel < 2) || (roi_lower <= roiLowerEnergy) );
    assert( (roi_upper >= roiUpperEnergy) || ((upper_channel + 3) >= dataH->num_gamma_channels()) );

    if( coFitPeaks.size() >= (upper_channel - lower_channel) ) //PeakFitDiffCostFunction will check for this in an assert
      throw runtime_error( "Invalid energy channels (" + std::to_string(lower_channel) + ", "
                          + std::to_string(upper_channel) + ") for " + std::to_string(coFitPeaks.size()) + " peaks." );
    assert( coFitPeaks.size() < (upper_channel - lower_channel) );

    const size_t nchannels = dataH->num_gamma_channels();
    const size_t midbin = dataH->find_gamma_channel( 0.5*(roiLowerEnergy + roiUpperEnergy) );
    const float binwidth = dataH->gamma_channel_width( midbin );
    const size_t num_fit_peaks = coFitPeaks.size() + 1;

    const PeakContinuum::OffsetType offset_type = input_continuum->type();
    const size_t num_continuum_pars = PeakContinuum::num_parameters( offset_type );

    const double prev_ref_energy = input_continuum->referenceEnergy();
    const bool prev_ref_energy_valid = ((prev_ref_energy >= roiLowerEnergy) && (prev_ref_energy <= roiUpperEnergy));
    const double reference_energy = prev_ref_energy_valid ? prev_ref_energy : roiLowerEnergy;

    const vector<bool> parFitFors = input_continuum->fitForParameter();
    assert( parFitFors.size() <= num_continuum_pars );

    const PeakDef::SkewType skew_type = PeakFitDiffCostFunction::skew_type_from_prev_peaks( coFitPeaks );
    const size_t num_skew = PeakDef::num_skew_parameters( skew_type );


    auto cost_functor = new PeakFitDiffCostFunction( dataH, coFitPeaks, roiLowerEnergy, roiUpperEnergy,
                                                    offset_type, reference_energy, skew_type, isHPGe, fit_options );

    //Choosing 8 paramaters to include in the `ceres::Jet<>` is 4 peaks in ROI, which covers most cases
    //  without introducing a ton of extra overhead.
    auto cost_function = new ceres::DynamicAutoDiffCostFunction<PeakFitDiffCostFunction,8>( cost_functor );

    const size_t num_fit_pars = cost_functor->number_parameters();

    cost_function->AddParameterBlock( static_cast<int>(num_fit_pars) );
    cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );

    const PeakFitDiffCostFunction::ProblemSetup prob_setup = cost_functor->get_problem_setup();

    assert( prob_setup.m_parameters.size() == num_fit_pars );
    assert( prob_setup.m_lower_bounds.size() == num_fit_pars );
    assert( prob_setup.m_upper_bounds.size() == num_fit_pars );

    vector<double> parameters = prob_setup.m_parameters;
    double * const pars = &parameters[0];


    ceres::Problem problem;

    // A brief look at a dataset of ~4k HPGe spectra with known truth-value peaks areas shows
    //  that using a loss function doesnt seem improve outcomes, when measured by sucessful
    //  peak fits, and by comparison of fit to truth peak areas.
    ceres::LossFunction *lossfcn = nullptr;

    problem.AddResidualBlock( cost_function, lossfcn, pars );


    if( !prob_setup.m_constant_parameters.empty() )
    {
      ceres::Manifold *subset_manifold = new ceres::SubsetManifold( static_cast<int>(num_fit_pars), prob_setup.m_constant_parameters );
      problem.SetManifold( pars, subset_manifold );
    }

    for( size_t i = 0; i < num_fit_pars; ++i )
    {
      if( prob_setup.m_lower_bounds[i].has_value() )
      {
        assert( *prob_setup.m_lower_bounds[i] <= parameters[i] );
        problem.SetParameterLowerBound(pars, static_cast<int>(i), *prob_setup.m_lower_bounds[i] );
      }

      if( prob_setup.m_upper_bounds[i].has_value() )
      {
        assert( *prob_setup.m_upper_bounds[i] >= parameters[i] );
        problem.SetParameterUpperBound(pars, static_cast<int>(i), *prob_setup.m_upper_bounds[i] );
      }
    }//for( size_t i = 0; i < num_fit_pars; ++i )


    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG
    options.use_nonmonotonic_steps = true;
    options.max_consecutive_nonmonotonic_steps = 10;

    //options.max_num_consecutive_invalid_steps = 10; //default 5

    // Initial trust region was very coursely optimized using a fit over ~4k HPGe spectra, using the median peak area
    //  error from truth value (whose value was 1.65*sqrt(area)); but also the average error was also about optimal
    //  (value 3.12*sqrt(area)).
    options.initial_trust_region_radius = 35*(cost_functor->number_parameters() - prob_setup.m_constant_parameters.size());
    //options.initial_trust_region_radius = 1e4; //
    options.max_trust_region_radius = 1e16;

    // Minimizer terminates when the trust region radius becomes smaller than this value.
    options.min_trust_region_radius = 1e-32;
    // Lower bound for the relative decrease before a step is accepted.
    options.min_relative_decrease = 1e-3; //pseudo optimized based on success rate of fitting peaks - but unknown effect on accuracy of fits.
    options.max_num_iterations = 50000;
#ifndef NDEBUG
    options.max_solver_time_in_seconds = 300.0;
#else
    options.max_solver_time_in_seconds = 30.0;
#endif

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    options.minimizer_progress_to_stdout = true;
    options.logging_type = ceres::PER_MINIMIZER_ITERATION;
#else
    options.minimizer_progress_to_stdout = false;
    options.logging_type = ceres::SILENT;
#endif

    /** Termination when `(new_cost - old_cost) < function_tolerance * old_cost`
     Default 1e-9;  Setting this to 1e-6 from 1e-9 is a ~80x speedup, and slgiht success fitting more peaks.
     Values from 1E-5 to 1E-9 dont seem to have a big effect on accuracy, but 1E-7 had best results for a HPGe dataset

     Value      AvrgAccuracyNSigma   MedianAccuracyNSigma   CPU-Time
     1e-5      3.25068                           1.66096                             465.296 s
     1e-6      3.20428                           1.69388                             581.317 s
     1e-7      3.16456                           1.6934                               1864.5 s
     1e-8      3.29325                           1.6936                               16833.9 s
     */
    options.function_tolerance = 1e-7;
    options.parameter_tolerance = 1e-11; //Default value is 1e-8.  Using 1e-11, so its usually the function tolerance that terminates things.
    options.num_threads = 1; //Probably wont have much/any effect

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    //std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to solve." << endl;
#endif
    
    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        break;
        
      case ceres::NO_CONVERGENCE:
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "The L-M ceres::Solver solving failed - NO_CONVERGENCE:\n" << summary.FullReport() << endl;
#endif
        throw runtime_error( "The L-M ceres::Solver solving failed - NO_CONVERGENCE." );

      case ceres::FAILURE:
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "The L-M ceres::Solver solving failed - FAILURE:\n" << summary.FullReport() << endl;
#endif
        throw runtime_error( "The L-M ceres::Solver solving failed - FAILURE." );
        
      case ceres::USER_FAILURE:
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "The L-M ceres::Solver solving failed - USER_FAILURE:\n" << summary.FullReport() << endl;
#endif
        throw runtime_error( "The L-M ceres::Solver solving failed - USER_FAILURE." );
    }//switch( summary.termination_type )

    cost_functor->m_ncalls = 0;

    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::DENSE_SVD; //SPARSE_QR;

    // Some terms we are fitting may not matter much, so when computing the inverse of J'J, we will opt
    //  to drop all terms where the eigen value of that term, divided by the maximum eigenvalue is
    //  less than the condition number.
    //  TODO: Currently leaving condition number at default 1e-14 - but what value is reasonable should be investigated
    cov_options.null_space_rank = -1;
    cov_options.min_reciprocal_condition_number = 1e-14;

    vector<double> uncertainties( num_fit_pars, 0.0 );
    double *uncertainties_ptr = uncertainties.data();
    
    
    ceres::Covariance covariance(cov_options);
    vector<pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back( make_pair( pars, pars) );
    
    if( !covariance.Compute(covariance_blocks, &problem) )
    {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
      cerr << "Failed to compute covariance!" << endl;
#endif
      uncertainties_ptr = nullptr;
    }else
    {
      vector<double> row_major_covariance( num_fit_pars * num_fit_pars );
      const vector<const double *> const_par_blocks( 1, pars );

      const bool success = covariance.GetCovarianceMatrix( const_par_blocks, row_major_covariance.data() );
      assert( success );
      if( success )
      {
        //Obviously we dont need to recopy things, just to get the diagonal, but eventually we will save the full matrix....
        vector<vector<double>> covariance_matrix( num_fit_pars, vector<double>(num_fit_pars, 0.0) );

        for( size_t row = 0; row < num_fit_pars; ++row )
        {
          for( size_t col = 0; col < num_fit_pars; ++col )
            covariance_matrix[row][col] = row_major_covariance[row*num_fit_pars + col];
        }//for( size_t row = 0; row < num_fit_pars; ++row )

        for( size_t i = 0; i < num_fit_pars; ++i )
        {
          // TODO: compare uncertainties with Minuit method - should maybe check out.
          if( covariance_matrix[i][i] > 0.0 )
            uncertainties[i] = sqrt( covariance_matrix[i][i] );
        }
      }else
      {
        uncertainties_ptr = nullptr;
      }//
    }//if( we failed to computer covariance ) / else

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    // Using numeric differentiation, should only be a single call to get covariance.
    cout << "Took " << cost_functor->m_ncalls.load() << " calls to get covariance." << endl;
#endif

    vector<double> residuals( cost_functor->number_residuals(), 0.0 );
    auto final_peaks = cost_functor->parametersToPeaks<PeakDef,double>( &parameters[0], uncertainties_ptr, residuals.data() );

    vector<shared_ptr<const PeakDef>> results( final_peaks.size() );
    for( size_t i = 0; i < final_peaks.size(); ++i )
      results[i] = make_shared<PeakDef>( final_peaks[i] );

    return results;
  }catch( std::exception &e )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    cout << "fit_peak_for_user_click_LM caught: " << e.what() << endl;
#endif
    //assert( 0 );
    throw;
  }//try / catch

  assert( 0 );
  throw runtime_error( "shouldnt have gotten here" );
  return {};
}//void fit_peaks_in_roi_LM(...)



void fit_peak_for_user_click_LM( PeakShrdVec &results,
                             const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                             const vector<shared_ptr<const PeakDef>> &coFitPeaksInput,
                             const double mean0, const double sigma0,
                             const double area0,
                             const float roiLowerEnergy,
                             const float roiUpperEnergy,
                             const bool isHPGe )
{
  vector<shared_ptr<const PeakDef>> coFitPeaks = coFitPeaksInput;
  
  const PeakDef::SkewType skew_type = PeakFitDiffCostFunction::skew_type_from_prev_peaks( coFitPeaks );
  const size_t num_skew = PeakDef::num_skew_parameters( skew_type );

  std::shared_ptr<PeakDef> candidatepeak = std::make_shared<PeakDef>(mean0, sigma0, area0);

  //The below should probably go off the number of bins in the ROI
  const size_t num_fit_peaks = coFitPeaksInput.size() + 1;
  PeakContinuum::OffsetType offset_type = (num_fit_peaks < (isHPGe ? 3 : 2)) ? PeakContinuum::Linear : PeakContinuum::Quadratic;

  if( coFitPeaks.size() )
    offset_type = std::max( offset_type, coFitPeaks[0]->continuum()->type() );


  // Make sure all initial peaks share a continuum
  if( coFitPeaks.size() )
  {
    shared_ptr<PeakContinuum> initial_continuum = make_shared<PeakContinuum>( *coFitPeaks[0]->continuum() );
    initial_continuum->setRange( roiLowerEnergy, roiUpperEnergy );
    candidatepeak->setContinuum( initial_continuum );
    for( size_t i = 0; i < coFitPeaks.size(); ++i )
    {
      auto newpeak = make_shared<PeakDef>( *coFitPeaks[i] );
      newpeak->setContinuum( initial_continuum );
      coFitPeaks[i] = newpeak;
    }

    if( initial_continuum->type() != offset_type )
    {
      const vector<bool> parFitFors = initial_continuum->fitForParameter();
      assert( parFitFors.size() <= PeakContinuum::num_parameters( initial_continuum->type() ) );

      bool any_continuum_par_fixed = false;
      for( size_t fit_for_index = 0; fit_for_index < parFitFors.size(); ++fit_for_index )
        any_continuum_par_fixed |= !parFitFors[fit_for_index];

      if( !any_continuum_par_fixed )
      {
        const double prev_ref_energy = initial_continuum->referenceEnergy();
        const bool prev_ref_energy_valid = ((prev_ref_energy >= roiLowerEnergy) && (prev_ref_energy <= roiUpperEnergy));
        const double reference_energy = prev_ref_energy_valid ? prev_ref_energy : roiLowerEnergy;

        initial_continuum->calc_linear_continuum_eqn( dataH, reference_energy, roiLowerEnergy, roiUpperEnergy, 2, 2 );
      }
      initial_continuum->setType( offset_type );
    }
  }else
  {
    shared_ptr<PeakContinuum> initial_continuum = candidatepeak->continuum();
    initial_continuum->setRange( roiLowerEnergy, roiUpperEnergy );
    const double reference_energy = roiLowerEnergy;
    initial_continuum->calc_linear_continuum_eqn( dataH, reference_energy, roiLowerEnergy, roiUpperEnergy, 2, 2 );
    initial_continuum->setType( offset_type );
  }//if( coFitPeaks.size() )

  coFitPeaks.push_back( candidatepeak );
  std::sort( coFitPeaks.begin(), coFitPeaks.end(), &PeakDef::lessThanByMeanShrdPtr );


  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options( 0 );

  try
  {
    results = fit_peaks_in_roi_LM( coFitPeaks, dataH, isHPGe, fit_options );
  }catch( std::exception &e )
  {
    results.clear();
  }
}//void fit_peak_for_user_click_LM(...)


void fit_peaks_LM( vector<shared_ptr<const PeakDef>> &results,
                  const vector<shared_ptr<const PeakDef>> input_peaks,
                  shared_ptr<const SpecUtils::Measurement> data,
                  const double stat_threshold,
                  const double hypothesis_threshold,
                  const bool is_refit,
                  const bool isHPGe ) throw()
{
  try
  {
    // Check all input peaks share a ROI
    for( size_t i = 1; i < input_peaks.size(); ++i )
    {
      if( input_peaks[i]->continuum() != input_peaks[0]->continuum() )
        throw runtime_error( "fit_peaks_LM: all input peaks must share ROI." );
    }

    results.clear();

    //We have to separate out non-gaussian peaks since they cant enter the
    //  fitting methods
    vector<shared_ptr<const PeakDef>> near_peaks, datadefined_peaks;

    shared_ptr<const SpecUtils::Measurement> ext_continuum;

    near_peaks.reserve( input_peaks.size() );

    for( const shared_ptr<const PeakDef> &p : input_peaks )
    {
      if( p->gausPeak() )
        near_peaks.push_back( p );
      else
        datadefined_peaks.push_back( p );

      if( !!p->continuum()->externalContinuum() )
        ext_continuum = p->continuum()->externalContinuum();
    }//for( const PeakDef &p : all_near_peaks )


    if( near_peaks.empty() )
      return;

    // Completely replace contents of `near_peaks` with new peaks, and new `PeakContinuum`
    //  so we dont affect the input peaks.
    local_unique_copy_continuum( near_peaks );

    //Need to make sure near_peaks and fixedpeaks are all gaussian (if not
    //  seperate them out, and add them in later).  If fitpeaks is non-gaussian
    //  ignore it or throw an exception.

    double lowx( 0.0 ), highx( 0.0 );

    {
      const shared_ptr<const PeakDef> &lowgaus = near_peaks.front();
      const shared_ptr<const PeakDef> &highgaus = near_peaks.back();

      double dummy = 0.0;
      findROIEnergyLimits( lowx, dummy, *lowgaus, data, isHPGe );
      findROIEnergyLimits( dummy, highx, *highgaus, data, isHPGe );
    }


    Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options( 0 );
    if( is_refit )
      fit_options |= PeakFitLM::PeakFitLMOptions::MediumRefinementOnly;

    results = fit_peaks_in_roi_LM( near_peaks, data, isHPGe, fit_options );

    // I'm pretty sure things have been sorted, but lets check
    assert( std::is_sorted(begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr ) );

    //
    //set_chi2_dof( data, results, 0, near_peaks.size() );


    for( size_t i = 1; i <= results.size(); ++i ) //Note weird convntion of index
    {
      const shared_ptr<const PeakDef> &peak = (results[i-1]);

      // Dont remove peaks whos amplitudes we arent fitting
      if( !peak->fitFor(PeakDef::GaussAmplitude) )
        continue;

      // Dont enforce a significance test if we are refitting the peak, and
      //  Sigma and Mean are fixed - the user probably knows what they
      //  are are doing.
      if( is_refit
         && !peak->fitFor(PeakDef::Sigma)
         && !peak->fitFor(PeakDef::Mean) )
      {
        continue;
      }

      const double num_sigma = (peak->amplitude() / peak->amplitudeUncert());
      const bool is_sig = (stat_threshold <= 0.0) || (num_sigma >= stat_threshold);

      const double dummy_stat_thresh = 0.0;
      const bool significant = chi2_significance_test( *peak, dummy_stat_thresh, hypothesis_threshold, {}, data );
      if( !is_sig || !significant )
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cerr) << "\tPeak at mean=" << peak->mean()
        << "is being discarded for not being significant"
        << "\n";
#endif
        results.erase( results.begin() + --i );
      }//if( !significant )
    }//for( size_t i = 1; i < fitpeaks.size(); ++i )

    bool removed_peak = false;
    for( size_t i = 1; i < results.size(); ++i ) //Note weird convention of index
    {
      const shared_ptr<const PeakDef> &this_peak = results[i-1];
      const shared_ptr<const PeakDef> &next_peak = results[i+1-1];

      // Dont remove peaks whos amplitudes we arent fitting
      if( !this_peak->fitFor(PeakDef::GaussAmplitude) )
        continue;

      const double min_sigma = min( this_peak->sigma(), next_peak->sigma() );
      const double mean_diff = next_peak->mean() - this_peak->mean();

      //In order to remove a gaussian, the peaks must both be within a sigma
      //  of eachother.  Note that this proccess doesnt care about the widths
      //  of the gaussians because we are assuming that the width of the gaussian
      //  should only be dependant on energy, so should only have one width of
      //  gaussian for a given energy
      if( (mean_diff/min_sigma) < 1.0 ) //XXX 1.0 chosen arbitrarily, and not checked
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cerr) << "Removing duplicate peak at x=" << this_peak->mean() << " sigma="
        << this_peak->sigma() << " in favor of mean=" << next_peak->mean()
        << " sigma=" << next_peak->sigma() << "\n";
#endif

        removed_peak = true;

        //Delete the peak with the worst chi2
        if( this_peak->chi2dof() > next_peak->chi2dof() )
          results.erase( begin(results) + i - 1 );
        else
          results.erase( begin(results) + i );

        i = i - 1; //incase we have multiple close peaks in a row
      }//if( (mean_diff/min_sigma) < 1.0 ) / else
    }//for( size_t i = 1; i < fitpeaks.size(); ++i )

    if( removed_peak )
      results = fit_peaks_in_roi_LM( results, data, isHPGe, fit_options );

    if( datadefined_peaks.size() )
    {
      results.insert( end(results), begin(datadefined_peaks), end(datadefined_peaks) );
      std::sort( begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr );
    }

    return;
  }catch( std::exception &e )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    cerr << "fit_peaks_LM: caught exception '" << e.what() << "'" << endl;
#endif
    results.clear();
    return;
  }//try / catch


  //We will only reach here if there was no exception, so since never expect
  //  this to actually happen, just assign the results to be same as the input
  assert( 0 );
  results = input_peaks;
}//vector<PeakDef> fit_peaks_LM(...);


vector<shared_ptr<const PeakDef>> fit_peaks_in_range_LM( const double x0, const double x1,
                                      const double ncausalitysigma,
                                      const double stat_threshold,
                                      const double hypothesis_threshold,
                                      const std::vector<std::shared_ptr<const PeakDef>> input_peaks,
                                      const std::shared_ptr<const SpecUtils::Measurement> data,
                                      const bool isRefit,
                                      const bool isHPGe )
{
  if( !data || (x1 < x0) )
    return input_peaks;

  vector<shared_ptr<const PeakDef>> all_peaks = input_peaks;
  std::sort( begin(all_peaks), end(all_peaks), &PeakDef::lessThanByMeanShrdPtr );

  const vector<shared_ptr<const PeakDef>> peaks_in_range = peaksInRange( x0, x1, ncausalitysigma, all_peaks );
  vector<shared_ptr<const PeakDef>> peaks_not_in_range = all_peaks;

  peaks_not_in_range.erase( std::remove_if( begin(peaks_not_in_range), end(peaks_not_in_range),
    [&peaks_in_range]( const shared_ptr<const PeakDef> &p ) -> bool {
      // Return false if the peak is in the energy range of interest
      return std::find(begin(peaks_in_range), end(peaks_in_range), p) != end(peaks_in_range);
    } ),
    end(peaks_not_in_range) );

  assert( (peaks_not_in_range.size() + peaks_in_range.size()) == input_peaks.size() );


  // The returned pointers all point to the original peaks.
  const vector<vector<shared_ptr<const PeakDef>>> seperated_peaks
                                          = causilyDisconnectedPeaks( ncausalitysigma, false, peaks_in_range );

  //Fit each of the ranges
  vector<vector<shared_ptr<const PeakDef>>> fit_peak_ranges( seperated_peaks.size() );
  SpecUtilsAsync::ThreadPool threadpool;
  for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  {
    threadpool.post( boost::bind( &fit_peaks_LM,
                                 boost::ref(fit_peak_ranges[peakn]),
                                 boost::cref(seperated_peaks[peakn]),
                                 data,
                                 stat_threshold,
                                 hypothesis_threshold,
                                 isRefit,
                                 isHPGe ) );
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  threadpool.join();

  vector<shared_ptr<const PeakDef>> results = peaks_not_in_range;

  //put the fit peaks back into 'input_peaks' so we can return all the peaks
  //  passed in, not just ones in X range of interest
  for( size_t peakn = 0; peakn < fit_peak_ranges.size(); ++peakn )
  {
    const vector<shared_ptr<const PeakDef>> &fit_peaks = fit_peak_ranges[peakn];
    results.insert( end(results), begin(fit_peaks), end(fit_peaks) );
  }//for( size_t peakn = 0; peakn < fit_peak_ranges.size(); ++peakn )

  std::sort( begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr );

  //Now make sure peaks from two previously causally disconnected regions
  //  didn't migrate towards each other, causing the regions to become
  //  causally connected now
  bool migration = false;
  for( size_t peakn = 1; peakn < fit_peak_ranges.size(); ++peakn )
  {
    if( fit_peak_ranges[peakn-1].empty() || fit_peak_ranges[peakn].empty() )
      continue;

    const shared_ptr<const PeakDef> &last_peak = fit_peak_ranges[peakn-1].back();
    const shared_ptr<const PeakDef> &this_peak = fit_peak_ranges[peakn][0];
    migration |= PeakDef::causilyConnected( *last_peak, *this_peak, ncausalitysigma, false );
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )

  if( migration )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    cerr << "fitPeaksInRange(...)\n\tWarning: Migration happened!" << endl;
#endif

    return fit_peaks_in_range_LM( x0, x1, ncausalitysigma,
                           stat_threshold, hypothesis_threshold,
                                 results, data, isRefit, isHPGe );
  }//if( migration )

  //  cout << "Fit took: " << timer.format() << endl;

  return results;
}//vector<shared_ptr<const PeakDef>> fit_peaks_in_range_LM(...)

std::vector<std::shared_ptr<const PeakDef>> refitPeaksThatShareROI_LM(
                                   const std::shared_ptr<const SpecUtils::Measurement> &data,
                                   const std::shared_ptr<const DetectorPeakResponse> &detector,
                                   const std::vector<std::shared_ptr<const PeakDef>> &inpeaks,
                                   const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options )
{
  vector<shared_ptr<const PeakDef>> answer;

  try
  {
    if( inpeaks.empty() )
      return answer;

    std::shared_ptr<const PeakContinuum> origCont = inpeaks[0]->continuum();

    for( const shared_ptr<const PeakDef> &p : inpeaks )
      if( origCont != p->continuum() )
        throw runtime_error( "refitPeaksThatShareROI_LM: all input peaks must share a ROI" );


    //const bool isHPGe = PeakFitUtils::is_high_res( spec );
    const auto resType = PeakFitUtils::coarse_resolution_from_peaks(inpeaks);
    const bool isHPGe = (resType == PeakFitUtils::CoarseResolutionType::High);

    answer = fit_peaks_in_roi_LM( inpeaks, data, isHPGe, fit_options );


    //now we need to go through and make sure the peaks we're adding are both
    //  significant, and improve the chi2/dof from before the fit.
    const double lx = origCont->lowerEnergy();
    const double ux = origCont->upperEnergy();
    const int lower_channel = static_cast<int>( data->find_gamma_channel( lx ) );
    const int upper_channel = static_cast<int>( data->find_gamma_channel( ux ) );
    const int nbin = (upper_channel > lower_channel) ? (upper_channel - lower_channel) : 1;
    const double prechi2Dof = chi2_for_region( inpeaks, data, lower_channel, upper_channel ) / nbin;
    const double postchi2Dof = chi2_for_region( answer, data, lower_channel, upper_channel ) / nbin;


    if( prechi2Dof < postchi2Dof )
    {
      answer.clear();

      if( true )
      {
        const double ncausalitysigma = 0.0;
        const double stat_threshold  = 0.0;
        const double hypothesis_threshold = 0.0;

        const bool isRefit = true;
        const vector<shared_ptr<const PeakDef>> refit_peaks
                     = fit_peaks_in_range_LM( lx, ux, ncausalitysigma, stat_threshold, hypothesis_threshold,
                                             inpeaks, data, isRefit, isHPGe );


        if( refit_peaks.size() == inpeaks.size() )
        {
          const double refit_chi2Dof = chi2_for_region( refit_peaks, data, lower_channel, upper_channel ) / nbin;
          cout << "refitPeaksThatShareROI_LM: refit_chi2Dof=" << refit_chi2Dof << ", prechi2Dof=" << prechi2Dof << ", postchi2Dof=" << postchi2Dof << endl;
          if( refit_chi2Dof <= prechi2Dof )
          {
            //cout << "Using re-fit peaks!" << endl;
            answer = refit_peaks;
          }
        }//if( output_peak.size() == inpeaks.size() )
      }

      if( answer.empty() )
        return answer;
    }//if( prechi2Dof >= 0.99*postchi2Dof )


    for( const shared_ptr<const PeakDef> &peak : answer )
    {
      const double mean = peak->mean();
      const double sigma = peak->sigma();
      const double data_counts = SpecUtils::gamma_integral( data, mean-0.5*sigma, mean+0.5*sigma );
      const double area = peak->gauss_integral( mean-0.5*sigma, mean+0.5*sigma );

      if( area < 5.0 || area < sqrt(data_counts) )
      {
        answer.clear();
        break;
      }//if( area < 5.0 || area < sqrt(data) )

      if( answer.size() == 1 )
      {
        static int ntimes = 0;
        if( ntimes++ < 3 )
          cerr << "refitPeaksThatShareROI_LM: Should perform some additional chi2 checks for the "
          "case answer.first.size() == 1" << endl;
      }//if( answer.first.size() == 1 )
    }//for( const PeakDefShrdPtr peak : answer.first )

  }catch( std::exception & )
  {
    answer.clear();
    cerr << "refitPeaksThatShareROI_LM: failed to find a better fit" << endl;
  }

  return answer;
}//refitPeaksThatShareROI_LM(...)

}//namespace PeakFitLM
