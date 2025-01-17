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
#include <thread>
#include <chrono>
#include <vector>
#include <memory>
#include <optional>
#include <iostream>

#include <Wt/WApplication>

#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#endif

#include "Eigen/Dense"
#include "ceres/ceres.h"
#include "ceres/loss_function.h"

#include "SandiaDecay/SandiaDecay.h"

#include "SpecUtils/SpecFile.h"
#include "SpecUtils/StringAlgo.h"
#include "SpecUtils/Filesystem.h"
#include "SpecUtils/D3SpectrumExport.h"

#include "SpecUtils/SpecUtilsAsync.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/MaterialDB.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/PhysicalUnits.h"
#include "InterSpec/ReactionGamma.h"
#include "InterSpec/UserPreferences.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/MassAttenuationTool.h"
#include "InterSpec/DetectorPeakResponse.h"
#include "InterSpec/GammaInteractionCalc.h"

#include "InterSpec/RelActCalc_imp.hpp"

using namespace std;

namespace Eigen {
  
  /** To use `Jets<>` in Eigen matrices, `ceres` defines
   `template <T,int N> struct NumTraits<ceres::Jet<T, N>>{...}`, however, it doesnt include
   `infinity()` or `quiet_NaN()`, which we need, so need to create an even more specialized version of the
   struct, with these functions.
   */
  template <int N>
  struct NumTraits<ceres::Jet<double, N>> {
    using Real = ceres::Jet<double, N>;
    using NonInteger = ceres::Jet<double, N>;
    using Nested = ceres::Jet<double, N>;
    using Literal = ceres::Jet<double, N>;
    
    static typename ceres::Jet<double, N> dummy_precision() { return ceres::Jet<double, N>(1e-12); }
    static inline Real epsilon() { return Real(std::numeric_limits<double>::epsilon()); }
    static inline Real infinity(){ return Real( std::numeric_limits<double>::infinity() ); }
    static inline Real quiet_NaN(){ return Real( std::numeric_limits<double>::quiet_NaN() ); }
    static inline int digits10() { return std::numeric_limits<double>::digits10; }
    static inline int max_digits10() { return std::numeric_limits<double>::max_digits10; }
    static inline Real highest() { return Real((std::numeric_limits<double>::max)()); }
    static inline Real lowest() { return Real(-(std::numeric_limits<double>::max)()); }
    
    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned,
      ReadCost = 1,
      AddCost = 1,
      // For Jet types, multiplication is more expensive than addition.
      MulCost = 3,
      HasFloatingPoint = 1,
      RequireInitialization = 1
    };
    
    template <bool Vectorized>
    struct Div {
      enum {
#if defined(EIGEN_VECTORIZE_AVX)
        AVX = true,
#else
        AVX = false,
#endif
        
        // Assuming that for Jets, division is as expensive as
        // multiplication.
        Cost = 3
      };
    };
  };//struct NumTraits<ceres::Jet<double, N>>
}//namespace Eigen

namespace
{
struct DoWorkOnDestruct
{
  std::function<void()> m_worker;
  DoWorkOnDestruct( std::function<void()> &&worker ) : m_worker( std::move(worker) ){}
  ~DoWorkOnDestruct(){ if(m_worker) m_worker(); }
};//struct DoWorkOnDestruct

/** We can either use a residual to force the normalization of the R.E. curve to 1.0 at the lowest
 * energy.  Or we can manually force the average measured relative efficiency to 1.0.
 * This does not apply to `RelEffEqnForm::FramPhysicalModel`.
 */
#define USE_RESIDUAL_TO_BREAK_DEGENERACY 0

  
template<typename T>
void fit_rel_eff_eqn_lls_imp( const RelActCalc::RelEffEqnForm fcn_form,
                              const size_t order,
                              const vector<double> &energies,
                              const vector<T> &data_values,
                              const vector<T> &data_uncertainties_orig,
                              vector<T> &fit_pars,
                              vector<vector<T>> *covariance )
{
  using namespace std;
  using namespace ceres;
  
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  const vector<T> &data_uncertainties = data_uncertainties_orig;
  
  assert( !data_values.empty() );
  assert( energies.size() == data_values.size() );
  assert( energies.size() == data_uncertainties.size() );
  
  if( data_values.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no data points provided." );
  
  const int poly_terms = static_cast<int>(order) + 1;
  const int num_peaks = static_cast<int>( data_values.size() );
  
  
  Eigen::MatrixX<T> A( num_peaks, poly_terms );
  Eigen::VectorX<T> b( num_peaks );
  
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const double energy = energies[row];
    const T measured_rel_eff = data_values[row];
    T uncertainty( data_uncertainties[row] );
    
    // Letting negative relative efficiencies through doesnt feel right, but I guess we'll do
    //  it to not mess up the L-M fitting of relative activities....
    //if( measured_rel_eff <= 0.0 )
    //  throw runtime_error( "fit_rel_eff_eqn_lls: Measured relative efficiency for energy "
    //                      + to_string(energy) + " is invalid ("
    //                      + to_string(measured_rel_eff) + ")" );
    
    //  But we'll put our foot down for negative or zero uncertainties.
    if( uncertainty <= 0.0 )
      throw runtime_error( "fit_rel_eff_eqn_lls: Uncertainty for energy " + to_string(energy)
                          + " is invalid." );
    
    switch( fcn_form )
    {
      case RelActCalc::RelEffEqnForm::LnX:
      {
        //y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
        b(row) = measured_rel_eff / uncertainty;
        break;
      }
        
      case RelActCalc::RelEffEqnForm::LnY:
      case RelActCalc::RelEffEqnForm::LnXLnY:
      case RelActCalc::RelEffEqnForm::FramEmpirical:
      {
        //LnY:           y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
        //LnXLnY:        y = exp (a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
        //FramEmpirical: y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
        //
        // We'll take the log of each side of the equation, and then solve for the parameters
        //
        // Note: when `f = ln(x)`, then `uncert(f) = uncert(x)/x`
        
        uncertainty = uncertainty / measured_rel_eff;
        
        // Note that we get the same answer (for a few problems I checked) if we use the following
        //  approximation to estimate uncertainty.
        //if( measured_rel_eff < 2.0*uncertainty )
        //  uncertainty = (2.0*uncertainty/measured_rel_eff) * fabs( log(0.75*measured_rel_eff) - log(1.25*measured_rel_eff) );
        //else
        //  uncertainty = 2.0*fabs( std::log(measured_rel_eff - 0.25*uncertainty)
        //                             - std::log(measured_rel_eff + 0.25*uncertainty) );
        
        b(row) = log(measured_rel_eff) / uncertainty;
        
        break;
      }
        
      case RelActCalc::RelEffEqnForm::FramPhysicalModel:
      {
        assert( 0 );
        throw runtime_error( "fit_rel_eff_eqn_lls: FramPhysicalModel not supported." );
      }
    }//switch( fcn_form )
    
    
    for( int col = 0; col < poly_terms; ++col )
    {
      switch( fcn_form )
      {
        case RelActCalc::RelEffEqnForm::LnX:
        case RelActCalc::RelEffEqnForm::LnXLnY:
          //y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
          // and
          //ln(y) = a + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ...
          A(row,col) = pow(log(energy), static_cast<T>(col)) / uncertainty;
          break;
          
        case RelActCalc::RelEffEqnForm::LnY:
          //ln(y) = a + b*x + c/x + d/x^2 + e/x^3 + ...
          if( col == 0 )
            A(row,col) = 1.0 / uncertainty;
          else if( col == 1 )
            A(row,col) = energy / uncertainty;
          else
            A(row,col) = pow(energy, 1.0 - col) / uncertainty;
          break;
          
        case RelActCalc::RelEffEqnForm::FramEmpirical:
          //ln(y) = a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3
          if( col == 0 )
            A(row,col) = 1.0 / uncertainty;
          else if( col == 1 )
            A(row,col) = (1.0 / (energy*energy)) / uncertainty;
          else
            A(row,col) = pow(log(energy), col - 1.0) / uncertainty;
          break;
          
        case RelActCalc::RelEffEqnForm::FramPhysicalModel:
          assert( 0 );
          throw runtime_error( "fit_rel_eff_eqn_lls: FramPhysicalModel not supported." );
      }//switch( fcn_form )
    }//for( int col = 0; col < poly_terms; ++col )
  }//for( int col = 0; col < poly_terms; ++col )
  
  // TODO: determine if HouseholderQr or BDC SVD is better/more-stable/faster/whatever
  //const Eigen::VectorXd solution = A.colPivHouseholderQr().solve(b);
  
  const Eigen::BDCSVD<Eigen::MatrixX<T>> bdc = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::VectorX<T> solution = bdc.solve(b);
  
  assert( solution.size() == (order + 1) );
  
  fit_pars.resize( solution.size() );
  for( size_t i = 0; i <= order; ++i )
    fit_pars[i] = solution(i);
  
  // Only compute covariance if it is wanted
  if( covariance )
  {
    // TODO: I'm sure Eigen::BDCSVD has the uncertainty matrix in it somewhere already computed, but
    //       for the moment (See pg 796 in Numerical Recipes for hint - probably has something to do
    //       with bdc.singularValues(), bdc.matrixV(), or bdc.matrixU()) we'll just be dumb and do
    //       extra (unstable?) work.
    const Eigen::MatrixX<T> A_transpose = A.transpose();
    const Eigen::MatrixX<T> alpha = Eigen::Product<Eigen::MatrixX<T>,Eigen::MatrixX<T>>( A_transpose, A ); //A_transpose * A;
    const Eigen::MatrixX<T> C = alpha.inverse();
    
    assert( C.rows() == solution.size() );
    assert( C.cols() == solution.size() );
    
    covariance->resize( solution.size() );
    
    for( size_t i = 0; i <= order; ++i )
    {
      vector<T> &row = (*covariance)[i];
      row.resize( solution.size() );
      for( size_t j = 0; j <= order; ++j )
        row[j] = C(i,j);
    }//for( loop over coefficients index )
  }//if( covariance )
}//fit_rel_eff_eqn_lls_imp(...)


template<typename T>
void fit_rel_eff_eqn_lls_imp( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<std::string> &isotopes,
                           const std::vector<T> &rel_acts,
                           const std::vector<RelActCalcManual::GenericPeakInfo> &peak_infos,
                           std::vector<T> &fit_pars,
                           std::vector<std::vector<T>> *covariance )
{
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  assert( !isotopes.empty() );
  if( isotopes.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no isotopes specified." );
  
  const int poly_terms = static_cast<int>(order) + 1;
  const int num_peaks = static_cast<int>( peak_infos.size() );
  
  vector<double> energies( num_peaks, 0.0 );
  vector<T> meas_rel_eff( num_peaks, T(0.0) ), meas_rel_eff_uncert( num_peaks, T(0.0) );
  
  bool unweighted_fit = false;
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const RelActCalcManual::GenericPeakInfo &peak = peak_infos[row];
    
    // A basic sanity check that the uncertainty in counts isnt garbage.
    if( peak.m_counts_uncert < 0.0 )
      throw runtime_error( "fit_rel_eff_eqn_lls: peak counts uncertainty can not be <0" );
    
    // Check there is a non-zero peak counts uncertainty; if its zero, we'll (arbitrarily) restrict
    //  to doing an un-weighted fit.  We could accept any non-zero peak.m_base_rel_eff_uncert
    //  and compute things just fine, but this would be highly suspect that the user has messed
    //  up filling out peak information, so we'll throw an exception.
    if( (peak.m_counts_uncert == 0.0) && (peak.m_base_rel_eff_uncert != -1.0) )
      throw runtime_error( "fit_rel_eff_eqn_lls: you must either provide a non-zero peak counts"
                          " uncertainty, or perform a unweighted fit" );
    
    const double energy = peak.m_energy;
    const double counts = peak.m_counts;
    const double counts_uncert = peak.m_counts_uncert;
    
    T raw_rel_counts( 0.0 );
    
    for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
    {
      const auto iso_pos = std::lower_bound( std::begin(isotopes), std::end(isotopes), line.m_isotope );
      assert( iso_pos != std::end(isotopes) );
      
      if( iso_pos == std::end(isotopes) )
        throw std::logic_error( "fit_rel_eff_eqn_lls: missing nuclide" );
      
      const size_t iso_index = static_cast<size_t>( iso_pos - std::begin(isotopes) );
      const T rel_act_value = rel_acts[iso_index];
      
      raw_rel_counts += line.m_yield * rel_act_value;
    }//for( const GenericLineInfo &line : peak.m_source_gammas )
    
    
    T measured_rel_eff = counts / raw_rel_counts;
    T measured_rel_eff_uncert = counts_uncert / raw_rel_counts;
    
    // We will clamp rel eff to zero or above ... this is a workaround since Eigens LM doesnt
    //  seem to allow constraining parameter ranges.
    if( (measured_rel_eff <= static_cast<double>(numeric_limits<float>::epsilon()))
       || isinf(measured_rel_eff)
       || isinf(measured_rel_eff) )
    {
      measured_rel_eff = T( 0.0 );
      if( peak.m_base_rel_eff_uncert > static_cast<double>(numeric_limits<float>::epsilon()) )
        measured_rel_eff_uncert = T( 0.0 );
      else
        measured_rel_eff_uncert = T( 1.0 );
    }
    
    if( peak.m_base_rel_eff_uncert == -1.0 )
    {
      if( row && !unweighted_fit )
        throw runtime_error( "fit_rel_eff_eqn_lls: for unweighted fit, all peaks must specify m_base_rel_eff_uncert == -1" );
      unweighted_fit = true;
      measured_rel_eff_uncert = T( 1.0 );
    }else
    {
      if( unweighted_fit )
        throw runtime_error( "fit_rel_eff_eqn_lls: for unweighted fit, all peaks must specify m_base_rel_eff_uncert == -1" );
      
      if( (peak.m_base_rel_eff_uncert < 0.0) || (peak.m_base_rel_eff_uncert > 1.0) )
        throw runtime_error( "fit_rel_eff_eqn_lls: m_base_rel_eff_uncert must be in range [0,1]" );
      
      if( peak.m_base_rel_eff_uncert > 0.0 )
      {
        // We should to be consistent with #ManualGenericRelActFunctor::eval in how we compute the
        // uncertainty
        
        const T add_uncert( counts * peak.m_base_rel_eff_uncert );
        measured_rel_eff_uncert = sqrt( pow(counts_uncert,2.0) + pow(add_uncert, 2.0) );
        measured_rel_eff_uncert /= raw_rel_counts;
      }//if( peak.m_base_rel_eff_uncert > 0.0 )
      
      // else keep as counts_uncert / raw_rel_counts
    }//if( do unweighted fit ) / else
    
    
    energies[row] = energy;
    meas_rel_eff[row] = measured_rel_eff;
    meas_rel_eff_uncert[row] = measured_rel_eff_uncert;
  }//for( int col = 0; col < poly_terms; ++col )
  
  
#if( !USE_RESIDUAL_TO_BREAK_DEGENERACY )
  
#ifdef _MSC_VER
#pragma message( "Double check how measured rel eff are being pinned to 1.0 - is there a better way?  Probably is!" )
#else
#warning "Double check how measured rel eff are being pinned to 1.0 - is there a better way?  Probably is!"
#endif
  
  const T sum_re = std::accumulate( begin(meas_rel_eff), end(meas_rel_eff), T(0.0) ); //Previous to 20250110, a value of 1,0 was used to initialize accumulate - not sure want that was, should probably check on this again
  const T average_re = sum_re / static_cast<double>( meas_rel_eff.size() );
  //const double first_re = meas_rel_eff[0];
  for( T &re : meas_rel_eff )
  {
    re /= average_re;
    //re -= average_re;        // This could go negative
    //re -= (first_re - 1.0);  // Seems to
    //re -= (meas_rel_eff[meas_rel_eff.size()-1] - 1.0);
  }
#endif
  
  
  fit_rel_eff_eqn_lls_imp<T>( fcn_form, order,
                          energies, meas_rel_eff, meas_rel_eff_uncert,
                          fit_pars, covariance );
}//fit_rel_eff_eqn_lls_imp(...)

/*
template<typename T>
void fit_rel_eff_eqn_lls_imp( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<SandiaDecayNucRelAct<T>> &nuclides,
                           const double base_rel_eff_uncert,
                           const std::vector<std::shared_ptr<const PeakDef>> &peak_infos,
                           vector<T> &fit_pars,
                           std::vector<std::vector<T>> *covariance )
{
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  assert( !nuclides.empty() );
  if( nuclides.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no nuclides specified." );
  
  assert( !peak_infos.empty() );
  if( peak_infos.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no peaks specified." );
  
  // We will map from the peaks mean, to the total number of gammas that contribute to that peak
  map<double,double> energy_gammas_map;
  
  // Get peak energies and widths (normally width this is just 'sigma', but for non-Gaussian peaks
  //  its 0.25 of the ROI)
  set<double> energies_seen; //a very poor check that there arent duplicate peaks
  vector< pair<double,double> > energy_widths, energy_obs_counts, energy_obs_counts_uncert;
  for( const auto &p : peak_infos )
  {
    // To use non-Gaussian peaks we would need to pass in the shared_ptr<const Measurement> data...
    //  maybe later if it ever matters
    if( !p->gausPeak() )
      throw runtime_error( "fit_rel_eff_eqn_lls: non-Gaussian peaks not supported yet" );
    
    // Note that in GammaInteractionCalc::ShieldingSourceChi2Fcn::observedPeakEnergyWidths
    //  we use the assigned nuclides gamma energy, as the energy - here we are using the peak mean.
    //  TODO: - revisit ether to use peak mean or its nuclide gamma as the energy - after implementing the rest of the manual RelAct calc stuff.
    const double energy = p->mean();
    const double sigma = p->gausPeak() ? p->sigma() : 0.25*p->roiWidth();
    const double amp = p->amplitude();
    //const double amp = p->gausPeak() ? p->amplitude() : p->areaFromData(data);
    const double ampUncert = p->amplitudeUncert();
    
    if( energies_seen.count(energy) )
      throw runtime_error( "fit_rel_eff_eqn_lls: multiple peaks with same energy - not allowed." );
    energies_seen.insert( energy );
    
    energy_widths.push_back( {energy, sigma} );
    energy_obs_counts.push_back( {energy, amp} );
    energy_obs_counts_uncert.push_back( {energy, ampUncert} );
  }//for( const PeakDef &peak : peaks )
  
  // JIC the peaks werent sorted, sort by just energies (although we did check no duplicate energies
  //  but we'll play it safe)
  auto sortByFirstOnly = []( const pair<double, double> &lhs, const pair<double, double> &rhs ){
    return lhs.first < rhs.first;
  };
  
  std::stable_sort( begin(energy_widths), end(energy_widths), sortByFirstOnly );
  std::stable_sort( begin(energy_obs_counts), end(energy_obs_counts), sortByFirstOnly );
  std::stable_sort( begin(energy_obs_counts_uncert), end(energy_obs_counts_uncert), sortByFirstOnly );
  
  
  // Now we will go through and get the amplitude of gammas we expect to contribute to a single peak
  //  (there may be multiple gammas from the same nuclide, as well as multiple nuclides that
  //   contribute to a single observable peaks).
  //  We will select a 'cluster' sigma of 1.5; this is what Activity/Shielding fit uses, but I dont
  //  think this value was derived by anything more than "that seems about right", and I havent
  //  run into an obvious case where this is not correct.
  const double photopeakClusterSigma = 1.5;
  set<const SandiaDecay::Nuclide *> nuclides_seen;
  for( const auto &n : nuclides )
  {
    if( nuclides_seen.count(n.nuclide) )
      throw runtime_error( "fit_rel_eff_eqn_lls: input nuclides must be unique" );
    
    nuclides_seen.insert( n.nuclide );
    
    SandiaDecay::NuclideMixture mixture;
    mixture.addNuclideByActivity(n.nuclide, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits);
    
    const double energyToCluster = -1;
    // TODO: we could account for decays during the measurement, but would need realTime here
    const bool accountForDecayDuringMeas = false;
    const double realTime = -1;
    GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities( energy_gammas_map,
                                                                          energy_widths, mixture, n.rel_activity, n.age,
                                                                          photopeakClusterSigma, energyToCluster,
                                                                          accountForDecayDuringMeas, realTime, nullptr, nullptr );
  }//for( const auto &n : nuclides )
  
  // Convert energy_gammas_map to a vector for convenience
  vector<pair<double,double>> energy_gammas;
  for( const auto &ec : energy_gammas_map )
    energy_gammas.push_back( ec );
  
  assert( energy_gammas.size() == peak_infos.size() );
  assert( energy_gammas.size() == energy_widths.size() );
  assert( energy_gammas.size() == energy_obs_counts.size() );
  assert( energy_gammas.size() == energy_obs_counts_uncert.size() );
  
  
  // Now put all this info onto a form so we can call into fit_rel_eff_eqn_lls(...); there is a
  //  commented out implementation of not having to do this, yet another, transformation of
  //  information.
  
  double max_pred_counts = 0.0;
  for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
    max_pred_counts = std::max(max_pred_counts, energy_gammas[peak_index].second);
  
  // Instead of keeping counts from each nuclide for each peak separate, we summed all nuclides
  //  together for each peak - so here we'll only use a single "Effective" isotope.
  const vector<string> isotopes{ "Effective" };
  const vector<double> rel_acts( 1, max_pred_counts );
  vector<GenericPeakInfo> generic_peak_infos;
  
  for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
  {
    GenericPeakInfo peak;
    peak.m_mean = peak_infos[peak_index]->mean();
    peak.m_energy = energy_gammas[peak_index].first;
    peak.m_fwhm = 2.35482*energy_widths[peak_index].second;
    peak.m_counts = energy_obs_counts[peak_index].second;
    peak.m_counts_uncert = (energy_obs_counts_uncert[peak_index].second > 0.0)
    ? energy_obs_counts_uncert[peak_index].second
    : sqrt(peak.m_counts);
    peak.m_base_rel_eff_uncert = base_rel_eff_uncert;
    
    const double yield = energy_gammas[peak_index].second / max_pred_counts;
    peak.m_source_gammas.emplace_back( yield, "Generic" );
    
    generic_peak_infos.push_back( peak );
  }//for( size_t row = 0; row < num_peaks; ++row )
  
  return fit_rel_eff_eqn_lls_imp<T>( fcn_form, order, isotopes, rel_acts, generic_peak_infos,
                             fit_pars, covariance );
}//fit_rel_eff_eqn_lls_imp(...)
*/
  

/** Setups the parameters for a shielding.
  
 Note that the number of parameters used is variable; when AN is being fit, a parameter will be used, otherwise not.
 A parameter for AD is always used.
 This was a result of the development path (specifically trouble getting auto-differentiation to work, and holding and
 manifolds to behave), and I think could be removed so number of parameter is fixed, but maybe we'll wait until
 we make another change (like being able to fix ratios of nuclides, or something) to change this
*/
void setup_physical_model_shield_par_manual( vector<int> &constant_parameters,
                                            double * const pars,
                                            size_t &index,
                                            const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt )
{
  if( !opt || (!opt->material && (opt->atomic_number == 0.0) && !opt->fit_atomic_number) )
    return;
  
  opt->check_valid();
  
  if( opt->material )
  {
    if( opt->fit_atomic_number )
      throw runtime_error( "You can not fit AN when defining a material" );
  }else if( opt->fit_atomic_number )
  {
    double lower_an = opt->lower_fit_atomic_number;
    double upper_an = opt->upper_fit_atomic_number;
    if( (lower_an == upper_an) && (lower_an == 0.0) )
    {
      lower_an = 1.0;
      upper_an = 98.0;
    }
    
    double an = opt->atomic_number;
    if( an < 1.0 )
    {
      assert( an == 0.0 );
      an = 0.5*(opt->lower_fit_atomic_number + opt->upper_fit_atomic_number);
    }
    
    if( (an < 1) || (an > 98) || (an < lower_an) || (an > upper_an) )
      throw runtime_error( "Self-atten AN is invalid" );
    
    pars[index] = an / RelActCalc::ns_an_ceres_mult;
    index += 1; //Add parameter for AN
    
    if( (lower_an < 1) || (upper_an > 98) || (upper_an <= lower_an) )
      throw runtime_error( "Self-atten AN limits is invalid" );
  }else
  {
    double an = opt->atomic_number;
    if( (an < 1) || (an > 98) )
      throw runtime_error( "Self-atten fixed AN is invalid" );
    
    // We wont actually add/use a parameter for fixed AN
  }//if( opt.material ) / else if( opt->fit_atomic_number )
  
  const double max_ad = RelActCalc::PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2;
  double ad = opt->areal_density / PhysicalUnits::g_per_cm2;
  double lower_ad = opt->lower_fit_areal_density / PhysicalUnits::g_per_cm2;
  double upper_ad = opt->upper_fit_areal_density / PhysicalUnits::g_per_cm2;
  
  if( (lower_ad == upper_ad) && (lower_ad == 0.0) )
  {
    lower_ad = 0.0;
    upper_ad = max_ad;
  }
  
  if( (ad == 0.0) && opt->fit_areal_density )
  {
    //ad = 0.5*(lower_ad + upper_ad); //Something like 250 would be way too much
    ad = 2.5; // We want something away from zero, because Ceres doesnt like zero values much - 2.5 is arbitrary
  }
  
  if( (ad < 0.0) || (ad > max_ad) )
    throw runtime_error( "Self-atten AD is invalid" );
  
  pars[index] = ad;
  
  
  if( opt->fit_areal_density )
  {
    // Check for limits
    if( (lower_ad < 0.0) || (upper_ad > max_ad) || (lower_ad >= upper_ad) )
      throw runtime_error( "Self-atten AD limits is invalid" );
  }else
  {
    constant_parameters.push_back( static_cast<int>(index) );
  }
  
  index += 1; //Add parameter for AD, always
}//void setup_physical_model_shield_par_manual( ceres::Problem... )


  
/** Functor for minimizing the relative activities; relative efficiency is fit for each set of
 activities via the #fit_rel_eff_eqn_lls function.
 */
struct ManualGenericRelActFunctor  /* : ROOT::Minuit2::FCNBase() */
{
  /** The form of relative efficiency equation to use. */
  //const RelActCalc::RelEffEqnForm m_eqn_form;
  
  /** The order of relative efficiency equation to use (not equation will have one more than this
   value coefficients)
   */
  //const int m_eqn_order;
  
  /** All isotope relative activities will be fit for using LM, with the relative eff curve forced
   to have a value near 1.0 for the first peak.
   */
  std::vector<string> m_isotopes;
  
  /** We will first normalize relative efficiency for each isotope to a flat line at y = 1.0.
   We will then use L-M to fit the multiples of these values that yield the best answer; this is
   to keep the values being fit for roughly around 1.0.
   
   i.e., This vector contains the relative activities for the relative efficiency line of y = 1.0
  (independent of energy).
   */
  vector<double> m_rel_act_norms;
  
  /** The input peak information.
   
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
   There will be one more residual than the number of peaks (with the very last residual being the
   difference of the relative efficiency, at the lowest energy, from 1.0).
#endif
   */
  //std::vector<RelActCalcManual::GenericPeakInfo> m_peak_infos;
  
  /** The input for the manual relative efficiency calculation. */
  RelActCalcManual::RelEffInput m_input;

  /** Warnings from setting up the problem - does not include problems evaluating things. */
  std::vector<std::string> m_setup_warnings;
  
  /** just for debug purposes, we'll keep track of how many times the eval function gets called. */
  mutable std::atomic<size_t> m_ncalls;
  
  /** Constructor for this functior.
   
   Will throw exception on error.
   */
  ManualGenericRelActFunctor( const RelActCalcManual::RelEffInput &input )
  : 
    m_isotopes{},
    m_rel_act_norms{},
  //m_eqn_form( input.eqn_form ),
  //m_eqn_order( static_cast<int>(input.eqn_order) ),
  //m_peak_infos( input.peaks ),
    m_input( input ),
    m_setup_warnings{},
    m_ncalls( 0 )
  {
    if( m_input.peaks.size() < 2 )
      throw runtime_error( "You must use at least two peaks." );
    
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      if( m_input.eqn_order != 0 )
        throw runtime_error( "ManualGenericRelActFunctor: equation order must be 0 for FramPhysicalModel." );
      if( !m_input.phys_model_detector )
        throw runtime_error( "ManualGenericRelActFunctor: detector must be specified for FramPhysicalModel." );

      try
      {
        if( m_input.phys_model_self_atten )
          m_input.phys_model_self_atten->check_valid();  

        for( const auto &opt : m_input.phys_model_external_attens )
        {
          if( !opt )
            throw runtime_error( "ManualGenericRelActFunctor: external attenuation may not be nullptr for FramPhysicalModel." );
          opt->check_valid();
        }
      }catch( const std::exception &e )
      {
        throw runtime_error( "ManualGenericRelActFunctor: attenuation input is invalid: " + std::string(e.what()) );
      }
    }else
    {
      // Apply some sanity checks to the eqn_order.  Realistically eqn_order should probably be
      //  between 3 and 6, but we'll allow an arbitrary amount of slop here.
      if( (m_input.eqn_order < 0) || (m_input.eqn_order >= 10) )
        throw runtime_error( "ManualGenericRelActFunctor: equation order must be at least 1 and less than 10." );

#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
      std::sort( begin(m_input.peaks), end(m_input.peaks),
                []( const RelActCalcManual::GenericPeakInfo &lhs, const RelActCalcManual::GenericPeakInfo &rhs ) -> bool {
        return lhs.m_energy < rhs.m_energy;
        } );
#endif
    }//if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else

    const size_t num_peaks = m_input.peaks.size();
    
    // We will check that we arent being passed in a ridiculous number of peak.
    const size_t max_allowed_peaks = 2000; // Arbitrarily chosen value
    if( num_peaks > max_allowed_peaks )
      throw std::runtime_error( "ManualGenericRelActFunctor: equation order must be at least 1 and less than 10." );
    
    bool unweighted_fit = false;
    
    // We'll go through and make a unique list of isotopes being fit for in m_isotopes, sorting them
    //  alphabetically (alphabetical ordering is arbitrary, other than for the order of isotope
    //  parameters to not be dependent on input ordering; i.e., for debugging purposes)
    //
    //  TODO: we should maybe order the nuclides by a niave guess as to the highest activity nuclide
    //        maybe just take the largest peak, or something...
    for( size_t peak_index = 0; peak_index < m_input.peaks.size(); ++peak_index )
    {
      const RelActCalcManual::GenericPeakInfo &peak = m_input.peaks[peak_index];
      
      if( peak.m_source_gammas.empty() )
        throw std::runtime_error( "ManualGenericRelActFunctor: Peak at " + std::to_string(peak.m_mean)
                                 + " keV has no source gammas defined." );
      
      if( peak.m_counts <= 0.0 )
        throw std::runtime_error( "ManualGenericRelActFunctor: peak counts must be >0.0." );
      
      if( (peak.m_counts_uncert < 0.0)
         || ((peak.m_counts_uncert == 0.0) && (peak.m_base_rel_eff_uncert != -1.0)) )
        throw std::runtime_error( "ManualGenericRelActFunctor: peak count uncertainty must be >0.0." );
      
      if( peak.m_base_rel_eff_uncert == -1.0 )
      {
        if( peak_index && !unweighted_fit )
          throw runtime_error( "ManualGenericRelActFunctor: RelActCalcManual::GenericPeakInfo::m_base_rel_eff_uncert must consistently either have a value [0,1], or -1" );
        
        unweighted_fit = true;
      }else if( unweighted_fit || (peak.m_base_rel_eff_uncert < 0.0) || (peak.m_base_rel_eff_uncert > 1.0) )
      {
        throw runtime_error( "ManualGenericRelActFunctor: RelActCalcManual::GenericPeakInfo::m_base_rel_eff_uncert must consistently either have a value [0,1], or -1" );
      }//
      
      set<string> nuclides_in_peak;
      
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      {
        const std::string &iso = line.m_isotope;
        
        // We probably wont usually have more than one line from a nuclide for a peak -- but we
        //  could so we'll just print out a warning, rather than throwing a hard error or something.
        if( nuclides_in_peak.count(iso) )
        {
          m_setup_warnings.push_back( "peak contained multiple lines from same nuclide." );
          std::cerr << "Warning: ManualGenericRelActFunctor: peak contained multiple lines from same nuclide.\n";
        }
        nuclides_in_peak.insert( iso );
        
        // lower_bound returns iterator to the that is not less than (i.e. greater or equal to)
        const auto pos = std::lower_bound( std::begin(m_isotopes), std::end(m_isotopes), iso );
        
        if( (pos == std::end(m_isotopes)) || ((*pos) != iso) )
          m_isotopes.insert( pos, iso );
        
        //  TODO: look at the various decay chains to determine better minimal yield value.
        //    Note that numeric_limits<float> is 1.19209e-07 on my system, which the BR for
        //    Pu that are important can be significantly smaller than this - so for the moment we'll
        //    essentially use zero.
        const double min_allowable_yield = std::numeric_limits<float>::min();
        if( line.m_yield <= min_allowable_yield )
          throw std::runtime_error( "ManualGenericRelActFunctor: yields must be greater than zero." );
      }//for( loop over nuclides that contribute to this peak )
    }//for( loop over fit peaks )
    
    const size_t num_isotopes = m_isotopes.size();
    
    // We should be guaranteed to have at least one isotope by here, but we'll throw in a check to
    //  avoid static analysis warnings, or jic, or something
    assert( num_isotopes );
    if( num_isotopes < 1 )
      throw std::runtime_error( "ManualGenericRelActFunctor: no isotopes specified." );
  
    // We'll first estimate the relative activities using a flat line
    vector<double> dummy_rel_act_norm_uncerts;
    const RelActCalc::RelEffEqnForm est_eqn_form = m_input.eqn_form;
    
    if( est_eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      RelActCalcManual::fit_act_to_phys_rel_eff( m_input, m_isotopes, m_input.peaks,
                        m_rel_act_norms, dummy_rel_act_norm_uncerts );

    }else
    {
      const vector<double> flat_rel_eff_coefs{ (est_eqn_form==RelActCalc::RelEffEqnForm::LnX ? 1.0 : 0.0), 0.0 };
      RelActCalcManual::fit_act_to_rel_eff( est_eqn_form, flat_rel_eff_coefs,
                       m_isotopes, m_input.peaks,
                       m_rel_act_norms, dummy_rel_act_norm_uncerts );
    }

    assert( m_isotopes.size() == m_rel_act_norms.size() );
    
    for( size_t i = 0; i < m_rel_act_norms.size(); ++i )
    {
      if( (m_rel_act_norms[i] < 1.0) || IsInf(m_rel_act_norms[i]) || IsNan(m_rel_act_norms[i]) )
      {
        m_setup_warnings.push_back( "The initial activity estimate for " + m_isotopes[i]
                                    + " was " + std::to_string(m_rel_act_norms[i])
                                    + ", so will use 1.0 instead.");
        m_rel_act_norms[i] = 1.0;
      }//if( m_rel_act_norms[i] < 1.0 )
    }//for( size_t i = 0; i < m_rel_act_norms.size(); ++i )
    
    if( num_isotopes > m_input.peaks.size() )
      throw std::runtime_error( "ManualGenericRelActFunctor: you must have at least as many peaks as"
                               " parameters you are fitting for." );
  }//ManualGenericRelActFunctor constructor
  
  size_t number_residuals() const
  {
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      return m_input.peaks.size();
    return m_input.peaks.size() + 1;
#else
    return m_input.peaks.size();
#endif
  }//size_t number_residuals() const
  
  size_t iso_index( const std::string &iso ) const
  {
    const auto iso_pos = std::lower_bound( std::begin(m_isotopes), std::end(m_isotopes), iso );
    assert( iso_pos != std::end(m_isotopes) );
    
    if( iso_pos == std::end(m_isotopes) )
      throw std::logic_error( "ManualGenericRelActFunctor: missing nuclide" );
    
    return static_cast<size_t>( iso_pos - std::begin(m_isotopes) );
  }
  
  template<typename T>
  T relative_activity( const std::string &iso, const vector<T> &x ) const
  {
    const size_t index = iso_index( iso );
    assert( index < x.size() );
    return m_rel_act_norms[index] * x[index];
  }
  
  
  template<typename T>
  struct PhysModelRelEqnDef
  {
    shared_ptr<const DetectorPeakResponse> det;
    
    std::optional<RelActCalc::PhysModelShield<T>> self_atten;
    std::vector<RelActCalc::PhysModelShield<T>> external_attens;
    
    std::optional<T> hoerl_b, hoerl_c;
  };//struct PhysModelRelEqnDef
  
  
  template<typename T>
  static PhysModelRelEqnDef<T> make_phys_eqn_input( const RelActCalcManual::RelEffInput &input,
                                         const std::vector<T> &eqn_coefficients )
  {
    // See `setup_physical_model_shield_par_manual(...)` for details of equation coefficients, but
    //  the short of it is there will only be a parameter for AN if material is nullptr, and
    //  atomic number is being fit for.
    // AD will always have a parameter dedicated to it.
    // if a RelActCalc::PhysicalModelShieldInput is nullptr, or material is null and
    //  opt->atomic_number == 0.0, and fit AN is false, then there will be no parameters
    //  dedicated to the PhysicalModelShieldInput.
    //  (this design is a bit vestigial, and we *could* change things to always have a
    //   consistent parameter definition, but maybe its worth waiting until we make another
    //   change to do this)
    
    assert( input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel );
    if( input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
      throw std::logic_error( "make_phys_eqn_input called for non-physical model." );
    
    PhysModelRelEqnDef<T> answer;
    
    answer.det = input.phys_model_detector;
    
    size_t shield_index = 0;
    if( input.phys_model_self_atten
       && (input.phys_model_self_atten->material
           || ((input.phys_model_self_atten->atomic_number >= 1.0) && (input.phys_model_self_atten->atomic_number <= 98.0))
           || input.phys_model_self_atten->fit_atomic_number) )
    {
      RelActCalc::PhysModelShield<T> self_atten;

      self_atten.material = input.phys_model_self_atten->material;
      if( !self_atten.material )
      {
        if( input.phys_model_self_atten->fit_atomic_number )
        {
          assert( (shield_index + 1) <= eqn_coefficients.size() );
          if( (shield_index + 1) > eqn_coefficients.size() )
            throw logic_error( "make_phys_eqn_input: not enough input coefficients (1)" );
          
          const T an = eqn_coefficients[shield_index] * RelActCalc::ns_an_ceres_mult;
          shield_index += 1;
          assert( an >= 1.0 && an <= 98.0 );
          self_atten.atomic_number = an;
        }else
        {
          self_atten.atomic_number = T( input.phys_model_self_atten->atomic_number );
        }
      }//if( !self_atten.material )
      
      assert( (shield_index + 1) <= eqn_coefficients.size() );
      if( (shield_index + 1) > eqn_coefficients.size() )
        throw logic_error( "make_phys_eqn_input: not enough input coefficients (2)" );
      
      assert( eqn_coefficients[shield_index] >= 0.0 && eqn_coefficients[shield_index] <= 500.0 );
      self_atten.areal_density = eqn_coefficients[shield_index] * PhysicalUnits::g_per_cm2;
      shield_index += 1;
      
      answer.self_atten = std::move(self_atten);
    }//m_options.phys_model_self_atten
      
      
    for( size_t i = 0; i < input.phys_model_external_attens.size(); ++i )
    {
      const auto &a = input.phys_model_external_attens[i];
      if( !a->material && ((a->atomic_number < 1.0) || (a->atomic_number > 98)))
        continue;
      
      RelActCalc::PhysModelShield<T> atten;
      atten.material = a->material;
      if( !atten.material )
      {
        if( a->fit_atomic_number )
        {
          assert( (shield_index + 1) <= eqn_coefficients.size() );
          if( (shield_index + 1) > eqn_coefficients.size() )
            throw logic_error( "make_phys_eqn_input: not enough input coefficients (3)" );
          
          const T an = eqn_coefficients[shield_index]  * RelActCalc::ns_an_ceres_mult;
          assert( (an >= 1.0) && (an <= 98.0) );
          atten.atomic_number = an;
          shield_index += 1;
        }else
        {
          atten.atomic_number = T( a->atomic_number );
        }
      }//if( !atten.material )
      
      assert( (shield_index + 1) <= eqn_coefficients.size() );
      if( (shield_index + 1) > eqn_coefficients.size() )
        throw logic_error( "make_phys_eqn_input: not enough input coefficients (4)" );
      
      assert( eqn_coefficients[shield_index] >= 0.0 && eqn_coefficients[shield_index] <= 500.0 );
      atten.areal_density = eqn_coefficients[shield_index] * PhysicalUnits::g_per_cm2;
      shield_index += 1;
      
      answer.external_attens.push_back( std::move(atten) );
    }//for( loop over input.phys_model_external_attens )
      
    if( input.phys_model_use_hoerl )
    {
      assert( (shield_index + 2) <= eqn_coefficients.size() );
      if( (shield_index + 2) > eqn_coefficients.size() )
        throw logic_error( "make_phys_eqn_input: not enough input coefficients (5)" );
      
      answer.hoerl_b = eqn_coefficients[shield_index];
      shield_index += 1;
      answer.hoerl_c = eqn_coefficients[shield_index];
      shield_index += 1;
    }//if( input.phys_model_use_hoerl )
    
    assert( shield_index == eqn_coefficients.size() );
    if( shield_index != eqn_coefficients.size() )
      throw logic_error( "make_phys_eqn_input: number of equation coefficients mismatch" );
    
    return answer;
  }//PhysModelRelEqnDef make_phys_eqn_input( const std::vector<double> &eqn_coefficients ) const
  
  
  template<typename T>
  std::function<T(double)> make_rel_eff_fcn( const std::vector<T> &eqn_coefficients ) const
  {
    const RelActCalc::RelEffEqnForm eqn_form = m_input.eqn_form;
    if( eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      // Note: will take copy of `eqn_coefficients` to avoid potential life-time issues, but we
      //       could maybe instead take a reference, for all our current use cases.
      return [eqn_form,eqn_coefficients]( double energy ){
        return eval_eqn_imp( energy, eqn_form, eqn_coefficients.data(), eqn_coefficients.size() );
      };
    }
      
    PhysModelRelEqnDef input = make_phys_eqn_input( m_input, eqn_coefficients );
    
    
    return [input]( double energy ){
      return eval_physical_model_eqn_imp<T>( energy, input.self_atten, input.external_attens,
                                            input.det.get(), input.hoerl_b, input.hoerl_c );
    };
  }//std::function<T(double)> rel_eff_fcn( const std::vector<T> &x ) const
  
  
  template<typename T>
  void eval_internal_lls_rel_eff( const std::vector<T> &x, T *residuals ) const
  {
    using namespace std;
    using namespace ceres;
    
    m_ncalls += 1;
    
    assert( residuals );
    assert( m_input.eqn_order >= 0 );
    assert( !m_input.use_ceres_to_fit_eqn );
    assert( x.size() == m_isotopes.size() );
    assert( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel );
    
    const size_t num_eqn_pars = static_cast<size_t>( m_input.eqn_order + 1 );
    const size_t num_isotopes = m_isotopes.size();
    assert( num_isotopes >= 1 );
    
    const T * const pars = x.data();
    vector<T> rel_activities( num_isotopes );
    for( size_t i = 0; i < num_isotopes; ++i )
      rel_activities[i] = this->relative_activity(m_isotopes[i], x);
    
    vector<T> eqn_coefficients;
    vector<vector<T>> eqn_cov;
    fit_rel_eff_eqn_lls_imp<T>( m_input.eqn_form, m_input.eqn_order,
                                m_isotopes,
                                rel_activities,
                                m_input.peaks,
                                eqn_coefficients, &eqn_cov );
    
    assert( eqn_coefficients.size() == (m_input.eqn_order + 1) );
    assert( eqn_cov.size() == (m_input.eqn_order + 1) );
    
    
    std::function<T(double)> rel_eff_curve = make_rel_eff_fcn<T>( eqn_coefficients );
    
    
    for( size_t index = 0; index < m_input.peaks.size(); ++index )
    {
      const RelActCalcManual::GenericPeakInfo &peak = m_input.peaks[index];
      
      const T curve_val = rel_eff_curve( peak.m_energy );
      
      T rel_src_counts( 0.0 );
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      {
        const T rel_activity = relative_activity( line.m_isotope, x );
        rel_src_counts += rel_activity * line.m_yield;
      }//for( const GenericLineInfo &line : peak.m_source_gammas )
      
      
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
      if( (index == 0) && (m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel) )
      {
        // Pin the first relative efficiency to 1.0; this is to prevent degeneracy of LM and LLS
        //  both essentially fitting for the overall normalization
        const double weight = 10.0; //Fairly arbitrary - not sure what to actually use - need something so the uncerts on activities dont get huge (which they will with no weight)
        residuals[m_peak_infos.size()] = weight*((peak.m_counts / rel_src_counts) - 1.0);
      }
#endif
      
      
      // TODO: should we fold in the uncertainty from the relative efficiency equation into the below as well - I dont think so since we are treating them as nuisance parameters?
      
      // Note: we want to compute:
      //  `(rel_efficiency - curve_val) / rel_eff_uncertainty`,
      //  But this has some divide by zero issues when rel act is zero, but the below gives the same
      //  residual, but avoiding these issues.
      if( peak.m_base_rel_eff_uncert == -1.0 )
      {
        // We are doing an unweighted fit
        
        // Avoid dividing by zero, so make sure rel_src_counts isnt really close to zero.
        // TODO: - need to evaluate using a different formulation where \c rel_src_counts being zero isnt a problem (I think its fine, but need to check before making the change - and make sure it wont effect how we use the covariances).
        if( ((rel_src_counts < 1.0E-8) && (rel_src_counts < (1.0E-6*peak.m_counts)))
           || (rel_src_counts < static_cast<double>(std::numeric_limits<float>::epsilon())) )
        {
          rel_src_counts = T(1.0E-6) * peak.m_counts;
        }
        
        residuals[index] = (peak.m_counts / rel_src_counts) - curve_val;
      }else if( peak.m_base_rel_eff_uncert == 0.0 )
      {
        // We are not using a m_base_rel_eff_uncert value
        const T pred_counts = curve_val * rel_src_counts;
        residuals[index] = (peak.m_counts - pred_counts) / peak.m_counts_uncert;
      }else
      {
        // We are using a m_base_rel_eff_uncert value
        assert( peak.m_base_rel_eff_uncert <= 1.0 );
        
        const T pred_counts = curve_val * rel_src_counts;
        // Note: for `add_uncert` below, we are using peak.m_counts, but it *could* also be
        //       reasonable to use `rel_src_counts` (which I was doing pre 20220720) or
        //       even `pred_counts`.  This is maybe worth revisiting.  Note that if you change
        //       how things are calculated here, you should also be consistent with
        //       #fit_rel_eff_eqn_lls.
        const double add_uncert = peak.m_counts * peak.m_base_rel_eff_uncert;
        const double uncert = sqrt( pow(peak.m_counts_uncert,2.0) + pow(add_uncert,2.0) );
        
        residuals[index] = (peak.m_counts - pred_counts) / uncert;
      }
    }//for( loop over energies to evaluate at )
  }//eval_internal_lls_rel_eff(...)


  template<typename T>
  void eval_internal_nl_rel_eff( const std::vector<T> &x, T *residuals ) const
  {
    using namespace std;
    using namespace ceres;
    
    m_ncalls += 1;
    
    assert( residuals );

    const size_t num_isos = num_isotopes();
  
    std::function<T (double)> rel_eff_fcn;
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      size_t par_num = num_isos;
      
      assert( m_input.phys_model_detector );
      const DetectorPeakResponse * const det = m_input.phys_model_detector.get();
      
      std::optional<RelActCalc::PhysModelShield<T>> self_atten;
      
      assert( !m_input.phys_model_self_atten
             || m_input.phys_model_self_atten->material
              || m_input.phys_model_self_atten->fit_atomic_number
              || ((m_input.phys_model_self_atten->atomic_number >= 1.0) && (m_input.phys_model_self_atten->atomic_number <= 98.0)) );
      
      if( m_input.phys_model_self_atten
         && (m_input.phys_model_self_atten->material
             || m_input.phys_model_self_atten->fit_atomic_number
             || ((m_input.phys_model_self_atten->atomic_number > 0.99) && (m_input.phys_model_self_atten->atomic_number < 98.001))) )
      {
        RelActCalc::PhysModelShield<T> att;
        att.material = m_input.phys_model_self_atten->material;
        
        if( att.material )
        {
          att.atomic_number = T(0.0);
        }else if( m_input.phys_model_self_atten->fit_atomic_number )
        {
          T an = x[par_num] * RelActCalc::ns_an_ceres_mult;
          att.atomic_number = fmin( fmax(an, 1.0), 98.0);
          par_num += 1;
        }else
        {
          att.atomic_number = T( fmin( fmax(m_input.phys_model_self_atten->atomic_number, 1.0), 98.0) );
        }
        
        att.areal_density = x[par_num] * PhysicalUnits::g_per_cm2;
        par_num += 1;
        
        self_atten = std::move(att);
      }//if( use internal attenuation shielding )
      
      vector<RelActCalc::PhysModelShield<T>> external_attens;
      for( const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &ext_atten : m_input.phys_model_external_attens )
      {
        assert( ext_atten );
        if( !ext_atten )
          continue;
        
        assert( ext_atten->material
                || ext_atten->fit_atomic_number
                || ((ext_atten->atomic_number >= 1.0) && (ext_atten->atomic_number <= 98.0)) );
        
        if( !ext_atten->material
           && !ext_atten->fit_atomic_number
           && ((ext_atten->atomic_number < 0.999) || (ext_atten->atomic_number > 98.001)) )
        {
          assert( 0 );
          continue;
        }
        
        RelActCalc::PhysModelShield<T> att;
        att.material = ext_atten->material;
        
        if( att.material )
        {
          att.atomic_number = T(0.0);
        }else if( ext_atten->fit_atomic_number )
        {
          T an = x[par_num] * RelActCalc::ns_an_ceres_mult;
          att.atomic_number = fmin( fmax(an, 1.0), 98.0);
          par_num += 1;
        }else
        {
          assert( (ext_atten->atomic_number >= 1.0) && (ext_atten->atomic_number <= 98.0) );
          att.atomic_number = T( ext_atten->atomic_number );
        }
        
        att.areal_density = x[par_num] * PhysicalUnits::g_per_cm2;
        par_num += 1;
        
        external_attens.push_back( std::move(att) );
      }//for( loop over external attenuators )
      
      
      std::optional<T> b, c;
      if( m_input.phys_model_use_hoerl )
      {
        b = x[par_num];
        par_num += 1;
        c = x[par_num];
        par_num += 1;
      }
      
      rel_eff_fcn = [self_atten,external_attens, det, b, c]( double energy ){
        return eval_physical_model_eqn_imp<T>( energy, self_atten, external_attens, det, b, c );
      };
      
      assert( par_num == x.size() );
    }else
    {
      assert( x.size() == (num_isos + m_input.eqn_order + 1) );
      rel_eff_fcn = [&x, this, num_isos]( double energy ){
        return eval_eqn_imp( energy, m_input.eqn_form, &(x[num_isos]), m_input.eqn_order + 1 );
      };
    }
    
    for( size_t index = 0; index < m_input.peaks.size(); ++index )
    {
      const RelActCalcManual::GenericPeakInfo &peak = m_input.peaks[index];
      
      T curve_val = rel_eff_fcn( peak.m_energy );
      
      T rel_src_counts{}; //scalar part of Jet will be default constructed to 0.0
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      {
        const T rel_activity = relative_activity( line.m_isotope, x );
        rel_src_counts += rel_activity * line.m_yield;
      }//for( const GenericLineInfo &line : peak.m_source_gammas )
      
      
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
      if( (index == 0) && (m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel) )
      {
        // Pin the first relative efficiency to 1.0; this is to prevent degeneracy of LM and LLS
        //  both essentially fitting for the overall normalization
        const double weight = 10.0; //Fairly arbitrary - not sure what to actually use - need something so the uncerts on activities dont get huge (which they will with no weight)
        residuals[m_peak_infos.size()] = weight*((peak.m_counts / rel_src_counts) - 1.0);
      }
#endif
      
      
      // TODO: should we fold in the uncertainty from the relative efficiency equation into the below as well - I dont think so since we are treating them as nuisance parameters?
      
      // Note: we want to compute:
      //  `(rel_efficiency - curve_val) / rel_eff_uncertainty`,
      //  But this has some divide by zero issues when rel act is zero, but the below gives the same
      //  residual, but avoiding these issues.
      if( peak.m_base_rel_eff_uncert == -1.0 )
      {
        // We are doing an unweighted fit
        
        // Avoid dividing by zero, so make sure rel_src_counts isnt really close to zero.
        // TODO: - need to evaluate using a different formulation where \c rel_src_counts being zero isnt a problem (I think its fine, but need to check before making the change - and make sure it wont effect how we use the covariances).
        if( ((rel_src_counts < 1.0E-8) && (rel_src_counts < (1.0E-6*peak.m_counts)))
           || (rel_src_counts < numeric_limits<T>::epsilon()) )
        {
          rel_src_counts = exp(T{}) * 1.0E-6 * peak.m_counts;
        }
        
        residuals[index] = (peak.m_counts / rel_src_counts) - curve_val;
      }else if( peak.m_base_rel_eff_uncert == 0.0 )
      {
        // We are not using a m_base_rel_eff_uncert value
        const T pred_counts = curve_val * rel_src_counts;
        residuals[index] = (peak.m_counts - pred_counts) / peak.m_counts_uncert;
      }else
      {
        // We are using a m_base_rel_eff_uncert value
        assert( peak.m_base_rel_eff_uncert <= 1.0 );
        
        const T pred_counts = curve_val * rel_src_counts;
        // Note: for `add_uncert` below, we are using peak.m_counts, but it *could* also be
        //       reasonable to use `rel_src_counts` (which I was doing pre 20220720) or
        //       even `pred_counts`.  This is maybe worth revisiting.  Note that if you change
        //       how things are calculated here, you should also be consistent with
        //       #fit_rel_eff_eqn_lls.
        const double add_uncert = peak.m_counts * peak.m_base_rel_eff_uncert;
        const double uncert = sqrt( pow(peak.m_counts_uncert,2.0) + pow(add_uncert,2.0) );
        
        residuals[index] = (peak.m_counts - pred_counts) / uncert;
      }
    }//for( loop over energies to evaluate at )
  }//eval_internal_nl_rel_eff(...)

/*
  template<typename T, int N>
  void eval( std::vector<ceres::Jet<T,N>> x, ceres::Jet<T,N> *residuals ) const
  {
    
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      const size_t num_isos = num_isotopes();
      // TODO: make sure scalar portion of Jet is valid range - need to use proper logic for parameter location...
      //We'll give the minimizer a little wiggle room to go negative for AD, but not to much.
      //  However, we'll make sure the AD is non-negative.
      //assert( x[num_isos + 1] > -1.0E-6 );
      //x[num_isos + 1].a = std::max( 0.0, x[num_isos + 1].a );
      //for( size_t i = 0; i < m_input.phys_model_external_attens.size(); ++i )
      //{
      //  assert( x[num_isos + 2 + 2*i + 1] > -1.0E-6 );
      //  x[num_isos + 2 + 2*i + 1].a = std::max( 0.0, x[num_isos + 2 + 2*i + 1].a );
      //}
    }//if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    
    if( !m_input.use_ceres_to_fit_eqn )
      throw std::logic_error( "ManualGenericRelActFunctor: relative efficiency equation should be fit using Ceres when using auto-differentiation." );

    eval_internal_nl_rel_eff( x, residuals );
  }
  */

  template<typename T>
  void eval( std::vector<T> x, T *residuals ) const
  {
    if( m_input.use_ceres_to_fit_eqn )
      eval_internal_nl_rel_eff<T>( x, residuals );
    else
      eval_internal_lls_rel_eff<T>( x, residuals );
  }
  
  
  virtual double operator()( const std::vector<double> &x ) const
  {
    vector<double> residuals( number_residuals(), 0.0 );
    try
    {
      eval( x, residuals.data() );
    }catch( std::exception &e )
    {
      cerr << "ManualGenericRelActFunctor::operator() caught: " << e.what() << endl;
      return std::numeric_limits<double>::max();
    }
    
    double chi2 = 0.0;
    for( size_t i = 0; i < m_input.peaks.size(); ++i )
      chi2 += residuals[i]*residuals[i];
    
    return chi2;
  }//operator() - for minuit
  
  
  // For Minuit2
  virtual double Up() const
  {
    return 1.0;
  }
  
  size_t num_isotopes() const
  {
    return m_isotopes.size();
  }

  size_t num_parameters() const
  {
    size_t num_pars = m_isotopes.size();
    if( !m_input.use_ceres_to_fit_eqn )
      return num_pars;
    
    if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      const auto &self_atten = m_input.phys_model_self_atten;
      if( self_atten
         && (self_atten->material
             || ((self_atten->atomic_number >= 0.999) && (self_atten->atomic_number <= 98.001))
             || self_atten->fit_atomic_number) )
      {
        if( self_atten->fit_atomic_number )
          num_pars += 1; //AN is only a parameter if being fit
        num_pars += 1; //AD is always a paramater
      }
      
      for( const auto &atten : m_input.phys_model_external_attens )
      {
        assert( atten );
        if( !atten )
          continue;
        
        assert( atten->material || ((atten->atomic_number >= 0.999) && (atten->atomic_number <= 98.001)) );
        if( atten->material
           || ((atten->atomic_number >= 0.999) && (atten->atomic_number <= 98.001))
           || atten->fit_atomic_number )
        {
          if( atten->fit_atomic_number )
            num_pars += 1;
          num_pars += 1;
        }
      }
      
      if( m_input.phys_model_use_hoerl )
        num_pars += 2;
    }else
    {
      num_pars += m_input.eqn_order + 1;
    }
    
    return num_pars;
  }
  
  // The return value indicates whether the computation of the
  // residuals and/or jacobians was successful or not.
  template<typename T>
  bool operator()(T const* const* parameters, T* residuals) const 
  {
    try
    {
      vector<T> pars( parameters[0], parameters[0] + num_parameters() );
     
      eval( pars, residuals );
    }catch( std::exception &e )
    {
      cerr << "ManualGenericRelActFunctor::operator() caught: " << e.what() << endl;
      return false;
    }
    
    return true;
  };//bool operator() - for Ceres

};//class ManualGenericRelActFunctor
}//namespace




namespace RelActCalcManual
{
GenericLineInfo::GenericLineInfo()
: m_yield( std::numeric_limits<double>::quiet_NaN() ),
  m_isotope( "InvalidIsotope" )
{
}

GenericLineInfo::GenericLineInfo( const double yield, const std::string &isotope )
 : m_yield( yield ),
   m_isotope( std::move(isotope) )
{
}


GenericPeakInfo::GenericPeakInfo()
 : m_energy( 0.0 ),
   m_mean( 0.0 ),
   m_fwhm( 0.0 ),
   m_counts( 0.0 ),
   m_counts_uncert( 0.0 ),
   m_base_rel_eff_uncert( 0.0 ),
   m_source_gammas{}
{
}


void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form, const size_t order,
                             const vector<double> &energies,
                             const vector<double> &data_values,
                             const vector<double> &data_uncertainties_orig,
                             vector<double> &fit_pars,
                             vector<vector<double>> *covariance )
{
  fit_rel_eff_eqn_lls_imp( fcn_form, order, energies, data_values,
                            data_uncertainties_orig, fit_pars, covariance );
}
  
  
void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                           const size_t order,
                           const std::vector<std::string> &isotopes,
                           const std::vector<double> &rel_acts,
                           const std::vector<GenericPeakInfo> &peak_infos,
                           std::vector<double> &fit_pars,
                           std::vector<std::vector<double>> *covariance )
{
  fit_rel_eff_eqn_lls_imp( fcn_form, order, isotopes, rel_acts, peak_infos, fit_pars, covariance );
}
  
  
vector<GenericPeakInfo> add_nuclides_to_peaks( const std::vector<GenericPeakInfo> &peaks,
                                              const std::vector<SandiaDecayNuc> &nuclides,
                                              const double real_time,
                                              const double cluster_sigma )
{
  vector<GenericPeakInfo> answer = peaks;
  
  vector< pair<double,double> > energy_widths, energy_obs_counts, energy_obs_counts_uncert;
  for( const auto &p : peaks )
    energy_widths.push_back( {p.m_energy, p.m_fwhm / 2.35482} );
  
  set<const SandiaDecay::Nuclide *> nuclides_seen;
  for( const auto &n : nuclides )
  {
    if( nuclides_seen.count(n.nuclide) )
      throw runtime_error( "add_nuclides_to_peaks: input nuclides must be unique" );
    nuclides_seen.insert( n.nuclide );
    
    SandiaDecay::NuclideMixture mixture;
    mixture.addNuclideByActivity(n.nuclide, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits);
    
    // We will map from the peaks mean, to the total number of gammas that contribute to that peak,
    //  for this nuclide
    map<double,double> energy_gammas_map;
        
    if( n.correct_for_decay_during_meas && (real_time <= 0) )
      throw runtime_error( "add_nuclides_to_peaks: measurement time must be specified if"
                          " correcting activities for nuclide decays during measurement.");
    
    const double activity = 1.0;
    const bool decay_correct = n.correct_for_decay_during_meas;
    
    GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities( energy_gammas_map,
                            energy_widths, mixture, activity, n.age, cluster_sigma, -1,
                            decay_correct, real_time, nullptr, nullptr );
    
    // Convert energy_gammas_map to a vector for convenience
    vector<pair<double,double>> energy_gammas;
    for( const auto &ec : energy_gammas_map )
      energy_gammas.push_back( ec );
    
    assert( energy_gammas.size() == answer.size() );
    
    for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
    {
      GenericPeakInfo &peak = answer[peak_index];
      const double yield = energy_gammas[peak_index].second;
      if( yield > numeric_limits<float>::min() )
        peak.m_source_gammas.emplace_back( yield, n.nuclide->symbol );
    }//for( size_t row = 0; row < num_peaks; ++row )
  }//for( const auto &n : nuclides )
  
  return answer;
}//add_nuclides_to_peaks(...)


/*
void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                         const size_t order,
                         const std::vector<SandiaDecayNucRelAct<double>> &nuclides,
                         const double base_rel_eff_uncert,
                         const std::vector<std::shared_ptr<const PeakDef>> &peak_infos,
                         vector<double> &fit_pars,
                         std::vector<std::vector<double>> *covariance )
{
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  assert( !nuclides.empty() );
  if( nuclides.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no nuclides specified." );
  
  assert( !peak_infos.empty() );
  if( peak_infos.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no peaks specified." );
  
  // We will map from the peaks mean, to the total number of gammas that contribute to that peak
  map<double,double> energy_gammas_map;
  
  // Get peak energies and widths (normally width this is just 'sigma', but for non-Gaussian peaks
  //  its 0.25 of the ROI)
  set<double> energies_seen; //a very poor check that there arent duplicate peaks
  vector< pair<double,double> > energy_widths, energy_obs_counts, energy_obs_counts_uncert;
  for( const auto &p : peak_infos )
  {
    // To use non-Gaussian peaks we would need to pass in the shared_ptr<const Measurement> data...
    //  maybe later if it ever matters
    if( !p->gausPeak() )
      throw runtime_error( "fit_rel_eff_eqn_lls: non-Gaussian peaks not supported yet" );
    
    // Note that in GammaInteractionCalc::ShieldingSourceChi2Fcn::observedPeakEnergyWidths
    //  we use the assigned nuclides gamma energy, as the energy - here we are using the peak mean.
    //  TODO: - revisit ether to use peak mean or its nuclide gamma as the energy - after implementing the rest of the manual RelAct calc stuff.
    const double energy = p->mean();
    const double sigma = p->gausPeak() ? p->sigma() : 0.25*p->roiWidth();
    const double amp = p->amplitude();
    //const double amp = p->gausPeak() ? p->amplitude() : p->areaFromData(data);
    const double ampUncert = p->amplitudeUncert();

    if( energies_seen.count(energy) )
      throw runtime_error( "fit_rel_eff_eqn_lls: multiple peaks with same energy - not allowed." );
    energies_seen.insert( energy );
    
    energy_widths.push_back( {energy, sigma} );
    energy_obs_counts.push_back( {energy, amp} );
    energy_obs_counts_uncert.push_back( {energy, ampUncert} );
  }//for( const PeakDef &peak : peaks )
  
  // JIC the peaks werent sorted, sort by just energies (although we did check no duplicate energies
  //  but we'll play it safe)
  auto sortByFirstOnly = []( const pair<double, double> &lhs, const pair<double, double> &rhs ){
    return lhs.first < rhs.first;
  };
  
  std::stable_sort( begin(energy_widths), end(energy_widths), sortByFirstOnly );
  std::stable_sort( begin(energy_obs_counts), end(energy_obs_counts), sortByFirstOnly );
  std::stable_sort( begin(energy_obs_counts_uncert), end(energy_obs_counts_uncert), sortByFirstOnly );
  
  
  // Now we will go through and get the amplitude of gammas we expect to contribute to a single peak
  //  (there may be multiple gammas from the same nuclide, as well as multiple nuclides that
  //   contribute to a single observable peaks).
  //  We will select a 'cluster' sigma of 1.5; this is what Activity/Shielding fit uses, but I dont
  //  think this value was derived by anything more than "that seems about right", and I havent
  //  run into an obvious case where this is not correct.
  const double photopeakClusterSigma = 1.5;
  set<const SandiaDecay::Nuclide *> nuclides_seen;
  for( const auto &n : nuclides )
  {
    if( nuclides_seen.count(n.nuclide) )
      throw runtime_error( "fit_rel_eff_eqn_lls: input nuclides must be unique" );
  
    nuclides_seen.insert( n.nuclide );
    
    SandiaDecay::NuclideMixture mixture;
    mixture.addNuclideByActivity(n.nuclide, GammaInteractionCalc::ShieldingSourceChi2Fcn::sm_activityUnits);
    
    const double energyToCluster = -1;
    // TODO: we could account for decays during the measurement, but would need realTime here
    const bool accountForDecayDuringMeas = false;
    const double realTime = -1;
    GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities( energy_gammas_map,
                                            energy_widths, mixture, n.rel_activity, n.age,
                                            photopeakClusterSigma, energyToCluster,
                                            accountForDecayDuringMeas, realTime, nullptr, nullptr );
  }//for( const auto &n : nuclides )
  
  // Convert energy_gammas_map to a vector for convenience
  vector<pair<double,double>> energy_gammas;
  for( const auto &ec : energy_gammas_map )
    energy_gammas.push_back( ec );

  assert( energy_gammas.size() == peak_infos.size() );
  assert( energy_gammas.size() == energy_widths.size() );
  assert( energy_gammas.size() == energy_obs_counts.size() );
  assert( energy_gammas.size() == energy_obs_counts_uncert.size() );
  
  
  // Now put all this info onto a form so we can call into fit_rel_eff_eqn_lls(...); there is a
  //  commented out implementation of not having to do this, yet another, transformation of
  //  information.
  
  double max_pred_counts = 0.0;
  for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
    max_pred_counts = std::max(max_pred_counts, energy_gammas[peak_index].second);
  
  // Instead of keeping counts from each nuclide for each peak separate, we summed all nuclides
  //  together for each peak - so here we'll only use a single "Effective" isotope.
  const vector<string> isotopes{ "Effective" };
  const vector<double> rel_acts( 1, max_pred_counts );
  vector<GenericPeakInfo> generic_peak_infos;
  
  for( size_t peak_index = 0; peak_index < energy_gammas.size(); ++peak_index )
  {
    GenericPeakInfo peak;
    peak.m_mean = peak_infos[peak_index]->mean();
    peak.m_energy = energy_gammas[peak_index].first;
    peak.m_fwhm = 2.35482*energy_widths[peak_index].second;
    peak.m_counts = energy_obs_counts[peak_index].second;
    peak.m_counts_uncert = (energy_obs_counts_uncert[peak_index].second > 0.0)
                             ? energy_obs_counts_uncert[peak_index].second
                             : sqrt(peak.m_counts);
    peak.m_base_rel_eff_uncert = base_rel_eff_uncert;
    
    const double yield = energy_gammas[peak_index].second / max_pred_counts;
    peak.m_source_gammas.emplace_back( yield, "Generic" );
    
    generic_peak_infos.push_back( peak );
  }//for( size_t row = 0; row < num_peaks; ++row )
  
  return fit_rel_eff_eqn_lls( fcn_form, order, isotopes, rel_acts, generic_peak_infos,
                              fit_pars, covariance );
}//fit_rel_eff_eqn_lls(...)
*/

void fit_act_to_rel_eff( const RelActCalc::RelEffEqnForm eqn_form,
                        const std::vector<double> &eqn_pars,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_uncerts )
{
  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    throw runtime_error( "fit_act_to_rel_eff: FramPhysicalModel not supported." );
  
  if( eqn_pars.empty() || eqn_pars.size() > 10 )
    throw runtime_error( "fit_act_to_rel_eff: invalid equation passed in." );

  auto eff_eqn = [eqn_form, eqn_pars]( double energy ){
    return eval_eqn( energy, eqn_form, eqn_pars );
  };
  
  fit_act_to_rel_eff( eff_eqn, isotopes, peak_infos, fit_rel_acts, fit_uncerts );
}//fit_act_to_rel_eff(...)


void fit_act_to_phys_rel_eff( const RelEffInput &input,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_uncerts )
{
  if( !input.phys_model_detector || !input.phys_model_detector->isValid() )
    throw runtime_error( "fit_act_to_phys_rel_eff: invalid detector passed in." );
  
  
  vector<int> dummy_const_pars;
  vector<double> rel_eff_pars( 2 + 2*input.phys_model_external_attens.size() + 2 ); //Likely too big - we'll resize later
  size_t rel_eff_index = 0;
  setup_physical_model_shield_par_manual( dummy_const_pars, rel_eff_pars.data(),
                                         rel_eff_index, input.phys_model_self_atten );
  
  for( const auto &opt : input.phys_model_external_attens )
  {
    assert( rel_eff_index < rel_eff_pars.size() );
    setup_physical_model_shield_par_manual( dummy_const_pars, rel_eff_pars.data(),
                                           rel_eff_index, opt );
  }
  
  if( input.phys_model_use_hoerl )
  {
    assert( rel_eff_index < rel_eff_pars.size() );
    rel_eff_pars[rel_eff_index] = 0.0; //(energy/1000)^b
    rel_eff_index += 1;
    assert( rel_eff_index < rel_eff_pars.size() );
    rel_eff_pars[rel_eff_index] = 1.0; //c^(1000/energy)
    rel_eff_index += 1;
  }//if( input.phys_model_use_hoerl )
  
  assert( rel_eff_index <= rel_eff_pars.size() );
  
  rel_eff_pars.resize( rel_eff_index );

  const ManualGenericRelActFunctor::PhysModelRelEqnDef eqn_input
        = ManualGenericRelActFunctor::make_phys_eqn_input( input, rel_eff_pars );
  
  const function<double(double)> eff_eqn
         = RelActCalc::physical_model_eff_function( eqn_input.self_atten, eqn_input.external_attens,
                                        eqn_input.det, eqn_input.hoerl_b, eqn_input.hoerl_c );

  fit_act_to_rel_eff( eff_eqn, isotopes, peak_infos, fit_rel_acts, fit_uncerts );
}//fit_act_to_rel_eff(...)


void fit_act_to_rel_eff( const std::function<double(double)> &eff_fcn,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_uncerts )
{
  // We want to fit the relative activities of the isotopes, given a relative efficiency curve
  
  // We'll start off with some basic sanity checks of the input.
  if( !eff_fcn )
    throw runtime_error( "fit_act_to_rel_eff: invalid rel eff fcn passed in." );
  
  if( peak_infos.size() < isotopes.size() )
    throw runtime_error( "fit_act_to_rel_eff: less peaks than isotopes." );
  
  std::vector<bool> used_isotopes( isotopes.size(), false );
  for( const GenericPeakInfo &info : peak_infos )
  {
    // We could blissfully ignore peaks with no source gammas (and assign them an activity of zero),
    //  but we'll be strict to maybe prevent some input mistakes, or something.
    if( info.m_source_gammas.empty() )
      throw runtime_error( "fit_act_to_rel_eff: peak at "
                          + to_string(info.m_mean) + " keV has no source gammas." );
    
    for( const GenericLineInfo &line : info.m_source_gammas )
    {
      if( line.m_yield <= 0.0 || IsInf(line.m_yield) || IsNan(line.m_yield) )
        throw runtime_error( "fit_act_to_rel_eff: invalid yield." );
      
      const auto pos = find( begin(isotopes), end(isotopes), line.m_isotope );
      if( pos == end(isotopes) )
        throw runtime_error( "fit_act_to_rel_eff: peak source isotope '" + line.m_isotope
                            + "' is not in list of wanted isotopes." );
      
      used_isotopes[pos - begin(isotopes)] = true;
    }//for( const GenericLineInfo &line : info.m_source_gammas )
    
    // Check peaks dont have the same nuclide multiple times
    for( size_t i = 1; i < info.m_source_gammas.size(); ++i )
    {
      for( size_t j = 0; j < i; ++j )
      {
        if( info.m_source_gammas[i].m_isotope == info.m_source_gammas[j].m_isotope )
          throw runtime_error( "fit_act_to_rel_eff: peak uses same isotope twice." );
      }//for( size_t j = 0; j < i; ++j )
    }//for( size_t i = 1; i < info.m_source_gammas.size(); ++i )
  }//for( const RelActCalc::GenericPeakInfo &info : peak_infos )
  
  assert( used_isotopes.size() == isotopes.size() );
  for( size_t i = 0; i < used_isotopes.size(); ++i )
  {
    if( !used_isotopes[i] )
      throw runtime_error( "fit_act_to_rel_eff: no peak with isotope '" + isotopes[i] + "' passed in." );
  }//for( size_t i = 0; i < used_isotopes.size(); ++i )
  
  // Checks are done, get to work
  
  //  We want to solve Ax = b, where
  //    Elements of A are the branching ratios for each isotope (e.g., where column correspond to
  //      isotopes), divided by the uncertainty of the peak - which corresponds to the row)
  //    x is the activities we are solving for
  //    b is the counts in each peak, divided by the relative efficiency for that energy (and all divided by the uncertainty of the peak).
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  
  const int num_isotopes = static_cast<int>( isotopes.size() );
  const int num_peaks = static_cast<int>( peak_infos.size() );
  
  Eigen::MatrixX<double> A = Eigen::MatrixX<double>::Zero( num_peaks, num_isotopes );
  Eigen::VectorX<double> b( num_peaks );
  
  
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const GenericPeakInfo &peak = peak_infos[row];
    const double energy = peak.m_energy;
    const double counts = peak.m_counts;
    const double counts_uncert = peak.m_counts_uncert;
    
    const double rel_eff = eff_fcn( energy );
    
    if( IsNan(rel_eff) || IsInf(rel_eff) )
      throw runtime_error( "fit_act_to_rel_eff: invalid rel eff at " + to_string(energy) + " keV." );
    
    const double rel_act = (counts / rel_eff);
    const double rel_act_uncert = rel_act * (counts_uncert / counts);
    
    b(row) = rel_act / rel_act_uncert;
    
    
    for( const GenericLineInfo &info : peak.m_source_gammas )
    {
      const auto pos = std::find( begin(isotopes), end(isotopes), info.m_isotope );
      assert( pos != end(isotopes) ); //we already checked this.
      
      const int column = static_cast<int>( pos - begin(isotopes) );
      assert( (column >= 0) && (column < static_cast<int>(isotopes.size())) ); //sometimes re-assurance is good
      
      // rel efficiency is just something like cps/rel_act
      assert( A(row,column) == 0.0 );
      A(row,column) = info.m_yield / rel_act_uncert;
    }//for( const GenericLineInfo &info : peak.m_source_gammas )
  }//for( size_t row = 0; row < num_peaks; ++row )
  
  // TODO: determine if HouseholderQr or BDC SVD is better/more-stable/faster/whatever
  //const Eigen::VectorXd solution = A.colPivHouseholderQr().solve(b);
  
  const Eigen::BDCSVD<Eigen::MatrixX<double>> bdc = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::VectorXd solution = bdc.solve(b);
  
  assert( solution.size() == num_isotopes );
  
  // TODO: I'm sure Eigen::BDCSVD has the uncertainty matrix in it somewhere already computed, but
  //       for the moment (See pg 796 in Numerical Recipes for hint - probably has something to do
  //       with bdc.singularValues(), bdc.matrixV(), or bdc.matrixU()) we'll just be dumb and do
  //       extra (unstable?) work.
  const Eigen::MatrixX<double> A_transpose = A.transpose();
  const Eigen::MatrixX<double> alpha = Eigen::Product<Eigen::MatrixX<double>,Eigen::MatrixX<double>>( A_transpose, A ); //A_transpose * A;
  const Eigen::MatrixX<double> C = alpha.inverse();
  
  fit_rel_acts.resize( solution.size() );
  fit_uncerts.resize( solution.size() );
  
  for( size_t i = 0; i < num_isotopes; ++i )
  {
    fit_rel_acts[i] = solution(i);
    fit_uncerts[i] = std::sqrt( C(i,i) );
  }
}//fit_act_to_rel_eff(...)


double RelEffSolution::relative_activity( const std::string &nuclide ) const
{
  for( const IsotopeRelativeActivity &i : m_rel_activities )
  {
    if( i.m_isotope == nuclide )
      return i.m_rel_activity;
  }
  
  throw runtime_error( "RelEffSolution::relative_activity: no nuclide '" + nuclide + "'" );
  return 0.0;
}//double relative_activity( const std::string &nuclide ) const


double RelEffSolution::relative_efficiency( const double energy ) const
{
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return eval_eqn( energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );

  const ManualGenericRelActFunctor::PhysModelRelEqnDef eqn_input
        = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );

  return RelActCalc::eval_physical_model_eqn( energy, eqn_input.self_atten, 
                              eqn_input.external_attens, eqn_input.det.get(), 
                              eqn_input.hoerl_b, eqn_input.hoerl_c );
}//double relative_efficiency( const double energy ) const


size_t RelEffSolution::nuclide_index( const std::string &nuc ) const
{
  for( size_t i = 0; i < m_rel_activities.size(); ++i )
  {
    if( m_rel_activities[i].m_isotope == nuc )
      return i;
  }
  
  throw runtime_error( "RelEffSolution::nuclide_index: invalid nuclide '" + nuc + "'" );
  return m_rel_activities.size(); //prevent warnings.
}//nuclide_index(...)


double RelEffSolution::activity_ratio( const size_t iso1, const size_t iso2 ) const
{
  assert( iso1 < m_rel_activities.size() );
  assert( iso2 < m_rel_activities.size() );
  return m_rel_activities[iso1].m_rel_activity / m_rel_activities[iso2].m_rel_activity;
}

double RelEffSolution::activity_ratio( const std::string &nuc1, const std::string &nuc2 ) const
{
  return activity_ratio( nuclide_index(nuc1), nuclide_index(nuc2) );
}


double RelEffSolution::activity_ratio_uncert( const size_t iso1, const size_t iso2 ) const
{
  assert( iso1 != iso2 );
  assert( iso1 < m_nonlin_covariance.size() );
  assert( iso2 < m_nonlin_covariance.size() );
  assert( iso1 < m_rel_activities.size() );
  assert( iso2 < m_rel_activities.size() );
  assert( m_nonlin_covariance.size() >= m_rel_activities.size() );
  assert( m_input.use_ceres_to_fit_eqn || (m_nonlin_covariance.size() == m_rel_activities.size()) );
  
  if( (iso1 == iso2)
     || (iso1 >= m_rel_activities.size())
     || (iso2 >= m_rel_activities.size()) )
  {
    //throw runtime_error( "RelEffSolution::activity_ratio_uncert: invalid iso number" );
    return -1;
  }
  
  const double norm_1 = m_activity_norms[iso1];
  const double norm_2 = m_activity_norms[iso2];
  
  // Divide the relative activities by their norms, to go back to the actual values fit for.
  const double fit_act_1 = m_rel_activities[iso1].m_rel_activity / norm_1;
  const double fit_act_2 = m_rel_activities[iso2].m_rel_activity / norm_2;
  
  const double cov_1_1 = m_nonlin_covariance[iso1][iso1];
  const double cov_1_2 = m_nonlin_covariance[iso1][iso2];
  const double cov_2_2 = m_nonlin_covariance[iso2][iso2];
  
  
  const double ratio = m_rel_activities[iso1].m_rel_activity
  / m_rel_activities[iso2].m_rel_activity;
  
  // TODO: I think this is the right way to compute ratio, taking into account correlations, given
  //       we actually fit for the values that multiplied m_activity_norms[],... need to double check
  const double uncert = fabs(ratio)
  * sqrt( (cov_1_1 / fit_act_1 / fit_act_1)
         + (cov_2_2 / fit_act_2 / fit_act_2)
         - (2.0 * cov_1_2 / fit_act_1 / fit_act_2) );
  
#ifndef NDEBUG
  {// Begin print out comparison between this result, and if we hadnt taken into account correlation
    const double iso1_uncert_frac = m_rel_activities[iso1].m_rel_activity_uncert / m_rel_activities[iso1].m_rel_activity;
    const double iso2_uncert_frac = m_rel_activities[iso2].m_rel_activity_uncert / m_rel_activities[iso2].m_rel_activity;
    const double niave_uncert = ratio * sqrt( iso1_uncert_frac*iso1_uncert_frac + iso2_uncert_frac*iso2_uncert_frac );
    
    cout
    << "Ratio " << m_rel_activities[iso1].m_isotope << "/" << m_rel_activities[iso2].m_isotope
    << " = " << ratio << " +- " << uncert << " (would be +- "
    << niave_uncert << " w/ no correlation)" << endl;
  }// End print out comparison between this result, and if we hadnt taken into account correlation
#endif
  
  return uncert;
}


double RelEffSolution::activity_ratio_uncert( const std::string &iso1, const std::string &iso2 ) const
{
  return activity_ratio_uncert( nuclide_index(iso1), nuclide_index(iso2) );
}


double RelEffSolution::mass_fraction( const std::string &nuclide ) const
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *wanted_nuc = db->nuclide( nuclide );
  if( !wanted_nuc )
    throw runtime_error( "RelEffSolution::mass_fraction('" + nuclide + "'): invalid nuclide" );
  
  double sum_rel_mass = 0.0, nuc_rel_mas = -1.0;
  for( size_t index = 0; index < m_rel_activities.size(); ++index )
  {
    const IsotopeRelativeActivity &act = m_rel_activities[index];
    const SandiaDecay::Nuclide * nuc = db->nuclide( act.m_isotope );
    if( !nuc )
      continue;
    
    const double rel_mass = act.m_rel_activity / nuc->activityPerGram();
    sum_rel_mass += rel_mass;
    
    if( nuc == wanted_nuc )
    {
      assert( nuc_rel_mas = -1.0 );
      nuc_rel_mas = rel_mass;
    }
  }//for( size_t index = 0; index < m_rel_activities.size(); ++index )
  
  if( nuc_rel_mas < 0.0 )
    throw runtime_error( "mass_fraction: invalid nuclide: " + nuclide );
  
  return nuc_rel_mas / sum_rel_mass;
}//double mass_fraction( const std::string &nuclide ) const

  
double RelEffSolution::mass_fraction( const std::string &nuclide, const double num_sigma ) const
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  const SandiaDecay::Nuclide *wanted_nuc = db->nuclide( nuclide );
  assert( wanted_nuc );
  if( !wanted_nuc )
    throw runtime_error( "RelEffSolution::mass_fraction('" + nuclide + "', num_sigma): invalid nuclide" );
  
  //assert( !m_nonlin_covariance.empty() ); // Failing to compute the activity covarances can happen sometimes
  if( m_nonlin_covariance.empty() )
    throw runtime_error( "RelEffSolution::mass_fraction('" + nuclide + "', num_sigma): no valid covariance." );
  
  assert( m_nonlin_covariance.size() >= m_activity_norms.size() );
  assert( m_input.use_ceres_to_fit_eqn || (m_nonlin_covariance.size() == m_rel_activities.size()) );
  
  const size_t nuc_index = std::find_if( begin(m_rel_activities), end(m_rel_activities),
                          [&nuclide]( const IsotopeRelativeActivity &val ) {
    return val.m_isotope == nuclide;
  }) - begin(m_rel_activities);
  
  assert( nuc_index < m_nonlin_covariance.size() );
  assert( nuc_index < m_rel_activities.size() );
  
  if( nuc_index >= m_rel_activities.size() )
    throw runtime_error( "RelEffSolution::mass_fraction('" + nuclide + "', "
                        + std::to_string(num_sigma) + "): nuclide not in solution set" );
  
  // The Covariance matrix is in terms of fit activities - however, we had divided out
  //  a normalization to bring them all near-ish 1.0 (before actually fitting for things
  //  though).
  const double norm_for_nuc = m_activity_norms[nuc_index];
  const double cov_nuc_nuc = m_nonlin_covariance[nuc_index][nuc_index];
  const double sqrt_cov_nuc_nuc = sqrt(cov_nuc_nuc);
  const double fit_act_for_nuc = m_rel_activities[nuc_index].m_rel_activity / norm_for_nuc;
  
  // Check that relative activity uncertainties have been computed compatible with what we are
  //  assuming here (and no funny business has happened).
  if( (norm_for_nuc * sqrt_cov_nuc_nuc) != m_rel_activities[nuc_index].m_rel_activity_uncert )
  {
    cout << "norm_for_nuc = " << norm_for_nuc << ", sqrt_cov_nuc_nuc = " << sqrt_cov_nuc_nuc << " (="<< norm_for_nuc*sqrt_cov_nuc_nuc << ")" << endl;
    cout << "m_rel_activities[nuc_index].m_rel_activity_uncert = " << m_rel_activities[nuc_index].m_rel_activity_uncert << endl;
    cout << "m_rel_activities[nuc_index].m_rel_activity = " << m_rel_activities[nuc_index].m_rel_activity << endl;
  }
  assert( (norm_for_nuc * sqrt_cov_nuc_nuc) == m_rel_activities[nuc_index].m_rel_activity_uncert );
  
  static std::atomic<int> num_times_here = 0;
  if( num_times_here < 5 )
  {
    ++num_times_here;
    cerr << "Do we need to modify how we compute mass-fraction uncertainties when `m_input.use_ceres_to_fit_eqn` is true?" << endl;
  }
  
  double sum_rel_mass = 0.0, nuc_rel_mas = -1.0;
  for( size_t index = 0; index < m_rel_activities.size(); ++index )
  {
    assert( m_nonlin_covariance[index].size() >= m_rel_activities.size() );
    
    const IsotopeRelativeActivity &act = m_rel_activities[index];
    const SandiaDecay::Nuclide * const nuc = db->nuclide( act.m_isotope );
    assert( nuc );
    if( !nuc )
      continue;
    
    const double norm_for_index = m_activity_norms[index];
    const double fit_act_for_index = m_rel_activities[index].m_rel_activity / norm_for_index;
    
    const double cov_nuc_index = m_nonlin_covariance[nuc_index][index];
    const double varied_fit_act_for_index = fit_act_for_index
                                  + (cov_nuc_index / cov_nuc_nuc) * num_sigma * sqrt_cov_nuc_nuc;
    
    const double rel_act = varied_fit_act_for_index * norm_for_index;
    const double rel_mass = rel_act / nuc->activityPerGram();
    
    sum_rel_mass += (std::max)( rel_mass, 0.0 );
    
    if( index == nuc_index )
    {
      assert( nuc == wanted_nuc );
      
      nuc_rel_mas = rel_mass;
    }//if( index == nuc_index )
  }//for( size_t index = 0; index < m_rel_activities.size(); ++index )
  
  if( nuc_rel_mas < 0.0 ) // This can happen when we go down a couple sigma
    return 0.0;
    
  return nuc_rel_mas / sum_rel_mass;
}//double mass_fraction( const std::string &iso, const double num_sigma ) const
    
  
std::string RelEffSolution::parameter_name( const size_t par_num ) const
{
  if( par_num >= m_fit_parameters.size() )
    throw runtime_error( "RelEffSolution::parameter_name: invalid parameter." );
  
  if( par_num < m_rel_activities.size() )
    return "Act(" + m_rel_activities[par_num].m_isotope + ")";
  
  if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    size_t working_par_num = m_rel_activities.size();
    if( working_par_num )
      working_par_num -= 1;
    
    if( m_input.phys_model_self_atten
       && (m_input.phys_model_self_atten->material
           || m_input.phys_model_self_atten->fit_atomic_number
           || ((m_input.phys_model_self_atten->atomic_number > 0.99) && (m_input.phys_model_self_atten->atomic_number < 98.001))) )
    {
      if( !m_input.phys_model_self_atten->material && m_input.phys_model_self_atten->fit_atomic_number )
      {
        working_par_num += 1;
        if( working_par_num == par_num )
          return "SAtt(AN)";
      }
      
      working_par_num += 1;
      if( working_par_num == par_num )
        return "SAtt(AD)";
    }//if( use internal attenuation shielding )
    
    for( size_t i = 0; i < m_input.phys_model_external_attens.size(); ++i )
    {
      const auto &ext_atten = m_input.phys_model_external_attens[i];
      if( !ext_atten )
        continue;
      
      if( !ext_atten->material
         && !ext_atten->fit_atomic_number
         && ((ext_atten->atomic_number < 0.999) || (ext_atten->atomic_number > 98.001)) )
      {
        continue;
      }
      
      if( !ext_atten->material && ext_atten->fit_atomic_number )
      {
        working_par_num += 1;
        if( working_par_num == par_num )
          return "EAtt" + std::to_string(i) + "(AN)";
      }
      
      working_par_num += 1;
      if( working_par_num == par_num )
        return "EAtt" + std::to_string(i) + "(AD)";
    }//for( loop over external attenuators )
    
    if( m_input.phys_model_use_hoerl )
    {
      working_par_num += 1;
      if( working_par_num == par_num )
        return "Hoerl(b)";
      working_par_num += 1;
      if( working_par_num == par_num )
        return "Hoerl(c)";
    }
    
    assert( 0 );
    throw std::logic_error( "Logic for determining Physical Model coefficient name is whack." );
  }//if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  
  
  return "RE_" + std::to_string( par_num - m_rel_activities.size() );
}//string parameter_name( const size_t par_num ) const
  

std::ostream &RelEffSolution::print_summary( std::ostream &strm ) const
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  switch( m_status )
  {
    case RelActCalcManual::ManualSolutionStatus::NotInitialized:
      strm << "Status: NotInitialized\n";
      break;
    case RelActCalcManual::ManualSolutionStatus::ErrorInitializing:
      strm << "Status: ErrorInitializing\n";
      break;
    case RelActCalcManual::ManualSolutionStatus::ErrorFindingSolution:
      strm << "Status: ErrorFindingSolution\n";
      break;
    case RelActCalcManual::ManualSolutionStatus::ErrorGettingSolution:
      strm << "Status: ErrorGettingSolution\n";
      break;
    case RelActCalcManual::ManualSolutionStatus::Success:
      strm << "Status: Success\n";
      break;
  }//switch( solution.m_status )
  
  if( !m_error_message.empty() )
  {
    strm << "--------------------------------------------------------------------------------\n";
    strm << "Error: " << m_error_message << "\n";
    strm << "--------------------------------------------------------------------------------\n";
  }
  
  if( !m_warnings.empty() )
  {
    strm << "--------------------------------------------------------------------------------\n";
    for( const string &warning : m_warnings )
      strm << "Warning: " << warning << "\n";
    strm << "--------------------------------------------------------------------------------\n";
  }//if( !m_warnings.empty() )
  
  strm << "Eqn coefficients: ";
  for( size_t i = 0; i < m_rel_eff_eqn_coefficients.size(); ++i )
  strm << (!i ? "" : ", ") << m_rel_eff_eqn_coefficients[i];
  strm << endl;
  strm << "Eqn coefficient covariance:\n";
  for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )
  {
    for( size_t j = 0; j < m_rel_eff_eqn_covariance[i].size(); ++j )
    strm << std::setw(14) << m_rel_eff_eqn_covariance[i][j];
    strm << endl;
  }
  strm << endl;
  
  strm << "Relative activities:" << endl;
  for( const RelActCalcManual::IsotopeRelativeActivity &i : m_rel_activities )
    strm << "\t" << i.m_isotope << ": " << i.m_rel_activity << " +- " << i.m_rel_activity_uncert << endl;
  
  
  double total_rel_mass = 0.0;
  for( const RelActCalcManual::IsotopeRelativeActivity &i : m_rel_activities )
  {
    const SandiaDecay::Nuclide *nuclide = db->nuclide( i.m_isotope );
    assert( nuclide );
    if( nuclide )
      total_rel_mass += i.m_rel_activity / nuclide->activityPerGram();
  }
  
  strm << "Relative masses:" << endl;
  for( const RelActCalcManual::IsotopeRelativeActivity &i : m_rel_activities )
  {
    const SandiaDecay::Nuclide *nuclide = db->nuclide( i.m_isotope );
    
    if( nuclide )
    {
      const double rel_mass = i.m_rel_activity / nuclide->activityPerGram();
      strm << "\t" << i.m_isotope << ": " << 100.0*rel_mass/total_rel_mass << endl;
    }
  }
  
  
  strm << "Chi2: " << m_chi2 << endl;
  strm << "Num Fcnt Evals: " << m_num_function_eval_total << "\n\n" << endl;
  
  
  for( const auto &peak: m_input.peaks )
  {
    double function_val = std::numeric_limits<double>::quiet_NaN();
    if( !m_rel_eff_eqn_coefficients.empty() )
    {
      if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        function_val = RelActCalc::eval_eqn( peak.m_energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );
      }else
      {
        const ManualGenericRelActFunctor::PhysModelRelEqnDef input
              = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );
        
        function_val = RelActCalc::eval_physical_model_eqn( peak.m_energy, input.self_atten,
                            input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_b );
      }//if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else
    }//if( !m_rel_eff_eqn_coefficients.empty() )

    strm << "For energy " << peak.m_energy << " (";
    for( size_t i = 0; i < peak.m_source_gammas.size(); ++i )
    strm << (i ? ", " : "") << peak.m_source_gammas[i].m_isotope;
    strm << ") function value is " << function_val << ", and:" << endl;
    
    //peak.m_counts
    strm << "\t";
    double total_contrib_counts = 0.0;
    for( size_t i = 0; i < peak.m_source_gammas.size(); ++i )
    {
      const RelActCalcManual::GenericLineInfo &line = peak.m_source_gammas[i];
      
      const string &iso = line.m_isotope;
      const double yield = line.m_yield;
      double rel_act = std::numeric_limits<double>::quiet_NaN();
      if( !m_rel_activities.empty() )
        rel_act = relative_activity(iso);
      
      const double contrib_counts = rel_act * yield * function_val;
      total_contrib_counts += contrib_counts;
      strm << (i ? ", " : "") << iso << ": fit " << contrib_counts << " counts";
    }//for( const RelActCalc::GammaLineInfo &line : peak.m_source_gammas )
    
    strm << " and observed " << peak.m_counts << "+-" << peak.m_counts_uncert << " (off by "
    << ((total_contrib_counts - peak.m_counts) / peak.m_counts_uncert)
    << " sigma)" << endl;
    
    
  }//for( const auto &peak: peak_infos )
  
  return strm;
}//std::ostream &print_summary( std::ostream &strm ) const


void RelEffSolution::get_mass_fraction_table( std::ostream &results_html ) const
{
  const int nsigfig = 4;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    return;
  
  results_html << "<table class=\"nuctable resulttable\">\n";
  results_html << "  <caption>Relative activities and mass fractions</caption>\n";
  results_html << "  <thead><tr>"
                 " <th scope=\"col\">Nuclide</th>"
                 " <th scope=\"col\">Rel. Act.</th>"
                 " <th scope=\"col\">Mass Frac.</th>"
                 //" <th scope=\"col\" title=\"This is percentage uncertainty - i.e., the percent error on the reported percent.\">Uncert.</th>"
                 " <th scope=\"col\" >2&sigma; range</th>"
                 " </tr></thead>\n";
  results_html << "  <tbody>\n";
  
  // We will put the normalized relative activity (i.e., activity divided by largest activity) in
  //  a tool-tip.
  //  However - we should probably just do this for the actual display - I'm just nervous of messing something up somewhere
  //  If we had more room, we could add a "Norm Act." column.
  double largest_rel_act = 0.0;
  for( const auto &act : m_rel_activities )
    largest_rel_act = std::max( largest_rel_act, act.m_rel_activity );
  
  for( size_t index = 0; index < m_rel_activities.size(); ++index )
  {
    const IsotopeRelativeActivity &act = m_rel_activities[index];
    const double uncert_percent = 100.0 * act.m_rel_activity_uncert / act.m_rel_activity;  //"percentage uncertainty"
    
    results_html << "  <tr><td>" << act.m_isotope << "</td>"
    << "<td title=\"The normalized relative activity (i.e., divided by largest rel. act.) is "
    << SpecUtils::printCompact( act.m_rel_activity / largest_rel_act, nsigfig + 1)
    << "\">" << SpecUtils::printCompact( act.m_rel_activity, nsigfig ) << "</td>";
    
    try
    {
      const double frac_mass = mass_fraction(act.m_isotope);
      results_html << "<td>" << SpecUtils::printCompact(100.0*frac_mass, nsigfig)      << "%</td>";
    }catch( std::exception & )
    {
      results_html << "<td>N.A.</td>";
    }
    
    // We need to do a better job of calculating mass fraction uncertainties.
    //  Maybe an okay way to go is increase activity by 1-sigma, then recalculate
    //  mass fraction, then re-normalize by new total mass
    
    string error_tt;
    try
    {
      const double frac_mass = mass_fraction(act.m_isotope);
      
      const double frac_mass_plus1 = mass_fraction(act.m_isotope, 1.0);
      const double frac_mass_minus1 = mass_fraction(act.m_isotope, -1.0);
       
      const double delta_plus1 = frac_mass_plus1 - frac_mass;
      const double delta_minus1 = frac_mass_minus1 - frac_mass;
      const double max_delta = (std::max)( fabs(delta_plus1), fabs(delta_minus1) );
      const double uncert_percent = 100.0 * max_delta / frac_mass;
      
      error_tt = "The 1-sigma mass fraction range is ["
      + SpecUtils::printCompact(100.0*frac_mass_minus1, nsigfig+1)
      + "% to "
      + SpecUtils::printCompact(100.0*frac_mass_plus1, nsigfig+1)
      + "%]. \nThe 1-sigma percentage uncertainty (i.e., the percent of the percent value) is "
      + SpecUtils::printCompact(uncert_percent, nsigfig+1)
      + "%";
    }catch( std::exception & )
    {
      // If we dont have covariance matrix, we will end up here
      error_tt = "No covariance available.";
    }
    
    results_html << "<td title=\"" << error_tt << "\">" ;
    //<< SpecUtils::printCompact(uncert_percent, nsigfig-1)
    try
    {
      const double frac_mass_plus2 = mass_fraction(act.m_isotope, 2.0);
      const double frac_mass_minus2 = mass_fraction(act.m_isotope, -2.0);
      
      const double frac_mass = mass_fraction(act.m_isotope);
      results_html << SpecUtils::printCompact(100.0*frac_mass_minus2, nsigfig-1)
      << "%, "
      << SpecUtils::printCompact( 100.0*frac_mass_plus2, nsigfig-1)
      << "%";
    }catch( std::exception & )
    {
      
    }
    
    results_html << "</td>"
    << "</tr>\n";
  }
  results_html << "  </tbody>\n"
  << "</table>\n\n";
}//void get_mass_fraction_table( std::ostream &strm ) const


void RelEffSolution::get_mass_ratio_table( std::ostream &results_html ) const
{
  const size_t nsigfig = 4;
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  if( !db )
    return;
  
  // Make the table of mass and activity ratios
  results_html << "<table class=\"massratiotable resulttable\">\n";
  results_html << "  <caption>Mass and Activity Ratios.</caption>\n";
  results_html << "  <thead><tr>"
  "<th scope=\"col\">Nuclides</th>"
  "<th scope=\"col\">Mass Ratio</th>"
  "<th scope=\"col\">Activity Ratio</th>"
  "</tr></thead>\n";
  results_html << "  <tbody>\n";
  
  for( size_t i = 1; i < m_rel_activities.size(); ++i )
  {
    for( size_t j = 0; j < i; ++j )
    {
      const string &nuc_i_str = m_rel_activities[i].m_isotope;
      const string &nuc_j_str = m_rel_activities[j].m_isotope;
      const SandiaDecay::Nuclide * const nuc_i = db->nuclide( nuc_i_str );
      const SandiaDecay::Nuclide * const nuc_j = db->nuclide( nuc_j_str );
      
      if( !nuc_i || !nuc_j )
      {
        results_html << "<tr><td>" << nuc_i_str << "/" << nuc_j_str << "</td><td>--</td><td>--</td></tr>\n";
        results_html << "<tr><td>" << nuc_j_str << "/" << nuc_i_str << "</td><td>--</td><td>--</td></tr>\n";
        continue;
      }
      
      const double act_i = relative_activity( nuc_i_str );
      const double act_j = relative_activity( nuc_j_str );
      
      const double mass_i = act_i / nuc_i->activityPerGram();
      const double mass_j = act_j / nuc_j->activityPerGram();
      
      const double i_to_j_specific_act = nuc_i->activityPerGram() / nuc_j->activityPerGram();
      const double j_to_i_specific_act = 1.0 / i_to_j_specific_act;
      
      const double i_to_j_act_ratio = activity_ratio( nuc_i_str, nuc_j_str );
      const double j_to_i_act_ratio = activity_ratio( nuc_j_str, nuc_i_str );
      
      const double i_to_j_mass_ratio = i_to_j_act_ratio * j_to_i_specific_act;
      const double j_to_i_mass_ratio = j_to_i_act_ratio * i_to_j_specific_act;
      
      if( !m_nonlin_covariance.empty() )
      {
        const double i_to_j_act_ratio_uncert = activity_ratio_uncert( nuc_i_str, nuc_j_str );
        const double j_to_i_act_ratio_uncert = activity_ratio_uncert( nuc_j_str, nuc_i_str );
        const double i_to_j_mass_ratio_uncert = i_to_j_act_ratio_uncert * j_to_i_specific_act;
        const double j_to_i_mass_ratio_uncert = j_to_i_act_ratio_uncert * i_to_j_specific_act;
        
        results_html << "<tr><td>" << nuc_i->symbol << "/" << nuc_j->symbol
        << "</td><td>" << PhysicalUnits::printValueWithUncertainty( i_to_j_mass_ratio, i_to_j_mass_ratio_uncert, nsigfig )
        << "</td><td>" << PhysicalUnits::printValueWithUncertainty( i_to_j_act_ratio, i_to_j_act_ratio_uncert, nsigfig )
        << "</td></tr>\n";
        
        results_html << "<tr><td>" << nuc_j->symbol << "/" << nuc_i->symbol
        << "</td><td>" << PhysicalUnits::printValueWithUncertainty( j_to_i_mass_ratio, j_to_i_mass_ratio_uncert, nsigfig )
        << "</td><td>"<< PhysicalUnits::printValueWithUncertainty( j_to_i_act_ratio, j_to_i_act_ratio_uncert, nsigfig )
        << "</td></tr>\n";
      }else
      {
        results_html << "<tr><td>" << nuc_i->symbol << "/" << nuc_j->symbol
        << "</td><td>" << SpecUtils::printCompact(i_to_j_mass_ratio, nsigfig) << " \xC2\xB1 n/a"
        << "</td><td>" << SpecUtils::printCompact(i_to_j_act_ratio, nsigfig) << " \xC2\xB1 n/a"
        << "</td></tr>\n";
        
        results_html << "<tr><td>" <<nuc_j->symbol << "/" << nuc_i->symbol
        << "</td><td>" << SpecUtils::printCompact(j_to_i_mass_ratio, nsigfig) << " \xC2\xB1 n/a"
        << "</td><td>" << SpecUtils::printCompact(j_to_i_act_ratio, nsigfig) << " \xC2\xB1 n/a"
        << "</td></tr>\n";
      }//if( we have covariance ) / else
    }//for( size_t j = 0; j < i; ++j )
  }//for( size_t i = 0; i < used_isotopes.size(); ++i )
  results_html << "  </tbody>\n"
  << "</table>\n\n";
}//void get_mass_ratio_table( std::ostream &strm ) const


double RelEffSolution::rel_eff_eqn_value( const double energy ) const
{
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return RelActCalc::eval_eqn( energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );
  
  const ManualGenericRelActFunctor::PhysModelRelEqnDef input 
                = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );

  return RelActCalc::eval_physical_model_eqn_imp<double>( energy, input.self_atten, input.external_attens,
                                            input.det.get(), input.hoerl_b, input.hoerl_c );
}


double RelEffSolution::rel_eff_eqn_uncert( const double energy ) const
{
  // Check if we use least linear squares to fit the equation; if so, use the covariance matrix directly
  if( !m_input.use_ceres_to_fit_eqn && (m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel) )
  {
    if( m_rel_eff_eqn_covariance.empty() )
      throw std::runtime_error( "RelEffSolution::rel_eff_eqn_uncert: Rel. Eff. Eqn. covariances not available." );

    const double val = RelActCalc::eval_eqn_uncertainty( energy, m_input.eqn_form,
                                            m_rel_eff_eqn_coefficients, m_rel_eff_eqn_covariance );
    if( isnan(val) )
      throw std::runtime_error( "RelEffSolution::rel_eff_eqn_uncert: NaN value for uncertainty." );

    return val;
  }//if( we can call eval_eqn_uncertainty )
  
  if( m_rel_eff_eqn_covariance.empty() )
    throw std::runtime_error( "RelEffSolution::rel_eff_eqn_uncert: nonlinear covariances not available." );

  assert( m_input.use_ceres_to_fit_eqn );

  assert( m_nonlin_covariance.size() == m_fit_parameters.size() );
  assert( m_rel_eff_eqn_covariance.size() == m_rel_eff_eqn_coefficients.size() );
  if( m_rel_eff_eqn_covariance.size() != m_rel_eff_eqn_coefficients.size() )
    throw std::logic_error( "RelEffSolution::rel_eff_eqn_uncert: covariance matrix does not match expected." );
  
  // I think we would be safe skipping this following check, at least on non-debug builds, but whatever
  for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )
  {
    assert( m_rel_eff_eqn_covariance[i].size() == m_rel_eff_eqn_covariance.size() );
    if( m_rel_eff_eqn_covariance[i].size() != m_rel_eff_eqn_covariance.size() )  //JIC for release builds
      throw runtime_error( "eval_eqn_uncertainty: covariance not a square matrix." );
  }//for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )

  double uncert_sq = 0.0;

  switch( m_input.eqn_form )
  {
    case RelActCalc::RelEffEqnForm::LnX:
    {
      // y = a + b*ln(x) + c*(ln(x))^2 + d*(ln(x))^3 + ...
      // We use/fit the same form where we use least linear squares, or Ceres, so we could just
      //  use the function we made for LLS fit
      const double log_energy = std::log(energy);
      for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )
      {
        for( size_t j = 0; j < m_rel_eff_eqn_covariance.size(); ++j )
          uncert_sq += std::pow(log_energy,1.0*i) * m_rel_eff_eqn_covariance[i][j] * std::pow(log_energy,1.0*j);
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
      assert( uncert_sq >= 0.0 );
      
      return sqrt(uncert_sq);
    }//case RelEffEqnForm::LnX:
      
      
    case RelActCalc::RelEffEqnForm::LnY:
    {
      // y = exp( a + b*x + c/x + d/x^2 + e/x^3 + ... )
      const double eval_val = eval_eqn( energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );
      
      const auto eval_derivative = [energy, eval_val]( const size_t order ) -> double {
        switch( order )
        {
          case 0:  return 1.0 * eval_val;
          case 1:  return energy * eval_val;
          default:
            break;
        }//switch( order )
        
        return eval_val * std::pow( energy, 1.0 - order );
      };//eval_derivative
      
      for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )
      {
        for( size_t j = 0; j < m_rel_eff_eqn_covariance.size(); ++j )
          uncert_sq += eval_derivative(i) * m_rel_eff_eqn_covariance[i][j] * eval_derivative(j);
      }//for( size_t i = 0; i < coefs.size(); ++i )
    }//case RelEffEqnForm::LnY:
      
    case RelActCalc::RelEffEqnForm::LnXLnY:
    {
      // y = exp(a  + b*(lnx) + c*(lnx)^2 + d*(lnx)^3 + ... )
      const double log_energy = std::log(energy);
      const double eval_val = eval_eqn( energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );
      
      for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )
      {
        for( size_t j = 0; j < m_rel_eff_eqn_covariance.size(); ++j )
          uncert_sq += (eval_val * std::pow(log_energy,1.0*i)) * m_rel_eff_eqn_covariance[i][j] * (eval_val * std::pow(log_energy,1.0*j));
      }//for( size_t i = 0; i < coefs.size(); ++i )
      
    }//case RelEffEqnForm::LnXLnY:
      
    case RelActCalc::RelEffEqnForm::FramEmpirical:
    {
      // y = exp( a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3 )
      const double log_energy = std::log(energy);
      const double eval_val = eval_eqn( energy, m_input.eqn_form, m_rel_eff_eqn_coefficients );
      
      double uncert_sq = 0.0;
      
      for( size_t i = 0; i < m_rel_eff_eqn_covariance.size(); ++i )
      {
        double i_component = 0.0;
        switch( i )
        {
          case 0:  i_component = eval_val;  break;
          case 1:  i_component = eval_val / (energy*energy); break;
          default: i_component = eval_val * std::pow(log_energy, i - 1.0); break;
        }//switch( i )
        
        for( size_t j = 0; j < m_rel_eff_eqn_covariance.size(); ++j )
        {
          double j_component = 0.0;
          switch( j )
          {
            case 0:  j_component = eval_val;  break;
            case 1:  j_component = eval_val / (energy*energy); break;
            default: j_component = eval_val * std::pow(log_energy, j - 1.0); break;
          }//switch( i )
          
          uncert_sq += i_component * m_rel_eff_eqn_covariance[i][j] * j_component;
        }
      }//for( size_t i = 0; i < coefs.size(); ++i )
    }//case RelEffEqnForm::FramEmpirical:
      
    case RelActCalc::RelEffEqnForm::FramPhysicalModel:
    {
      // y = [(1 - exp(-mu(AN_0,E)*x_j))/(mu(AN_0,E)*x_j)]
      //     * [exp(-mu(AN_1,E)*x_1) * exp(-mu(AN_2,E)*x_2) * ...]
      //     * [(E/1000)^b * c^(1000/E)]
      //     * [Det Eff]
      //
      // https://www.derivative-calculator.net
      // let f   = [(1 - exp(-function_m(x_0)*x_1))/(function_m(x_0)*x_1)]
      //           * [exp(-function_m(x_2)*x_3)*exp(-function_m(x_4)*x_5)]
      //           * ((E/1000)^(x_6)) * ((x_7)^(1000/E))
      //----------------------------------------------------------------------
      // df/dx_0 = ((E/1000)^x_6 * x_7^(1000 / E) * e^(-x_1 * function_m(x_0) - function_m(x_4) * x_5 - function_m(x_2) * x_3) * (x_1 * function_m(x_0) - e^(x_1 * function_m(x_0)) + 1) * function_m'(x_0))
      //    / (x_1 * function_m(x_0)^2)
      //
      // df/dx_1 = -(E^x_6 * x_7^(1000 / E) * (e^(function_m(x_0) * x_1) - function_m(x_0) * x_1 - 1) * e^(-function_m(x_0) * x_1 - function_m(x_4) * x_5 - function_m(x_2) * x_3))
      //    / (function_m(x_0) * 1000^x_6 * x_1^2)
      //
      // df/dx_2 = -(E^x_6 * (1 - e^(-function_m(x_0) * x_1)) * x_3 * x_7^(1000 / E) * e^(-x_3 * function_m(x_2) - function_m(x_4) * x_5) * function_m'(x_2))
      //    / (function_m(x_0) * x_1 * 1000^x_6)
      //
      // df/dx_3 = -(E^x_6 * (1 - e^(-function_m(x_0) * x_1)) * function_m(x_2) * x_7^(1000 / E) * e^(-function_m(x_2) * x_3 - function_m(x_4) * x_5))
      //    / (function_m(x_0) * x_1 * 1000^x_6)
      //
      // df/dx_4 = -(E^x_6 * (1 - e^(-function_m(x_0) * x_1)) * x_5 * x_7^(1000 / E) * e^(-x_5 * function_m(x_4) - function_m(x_2) * x_3) * function_m'(x_4))
      //    / (function_m(x_0) * x_1 * 1000^x_6)
      //
      // df/dx_5 = -(E^x_6 * (1 - e^(-function_m(x_0) * x_1)) * function_m(x_4) * x_7^(1000 / E) * e^(-function_m(x_4) * x_5 - function_m(x_2) * x_3))
      //    / (function_m(x_0) * x_1 * 1000^x_6)
      //
      // df/dx_6 = (E^x_6 * (ln(E) - ln(1000)) * (e^(function_m(x_0) * x_1) - 1) * e^(-function_m(x_4) * x_5 - function_m(x_2) * x_3 - function_m(x_0) * x_1) * x_7^(1000 / E))
      //    / (function_m(x_0) * x_1 * 1000^x_6)
      //
      // df/dx_7 = (E^(x_6 - 1) * (1 - e^(-function_m(x_0) * x_1)) * e^(-function_m(x_4) * x_5 - function_m(x_2) * x_3) * 1000^(1 - x_6) * x_7^(1000 / E - 1))
      //    / (function_m(x_0) * x_1)
      
      const auto is_valid_shield = []( const PhysModelShieldFit * const result,
                              const RelActCalc::PhysicalModelShieldInput * const input ) -> bool {
        if( !result || !input )
          return false;
        assert( (!result->m_material) == (!input->material) );
        if( result->m_material )
          return true;
        return (input->fit_atomic_number || ((result->m_atomic_number >= 1.0) && (result->m_atomic_number <= 98.0)));
      };//is_valid_shield(...)
      
      const auto an_was_fit = []( const RelActCalc::PhysicalModelShieldInput * const input ) -> bool {
        return input && input->fit_atomic_number;
      };
      
      const auto ad_was_fit = []( const RelActCalc::PhysicalModelShieldInput * const input ) -> bool {
        return input && input->fit_areal_density;
      };
      
      assert( (!m_phys_model_self_atten_shield) == (!m_input.phys_model_self_atten) );
      if( (!m_phys_model_self_atten_shield) != (!m_input.phys_model_self_atten) )
        throw runtime_error( "eval_eqn_uncertainty: result/input mismatch" );
      
      assert( m_phys_model_external_atten_shields.size() == m_input.phys_model_external_attens.size() );
      if( m_phys_model_external_atten_shields.size() != m_input.phys_model_external_attens.size() )
        throw runtime_error( "eval_eqn_uncertainty: result/input mismatch" );
      
      
      
      // Since the meaning of parameters doesnt map to their index, its a little awkward, we'll
      //  expand things so we can do a fixed mapping, but also track which variables we should skip.
      const size_t max_num_par = 2 + 2*m_phys_model_external_atten_shields.size() + 2;
      const size_t num_actual_pars = m_rel_eff_eqn_coefficients.size();
      vector<size_t> par_index( max_num_par, num_actual_pars ); //m
      vector<double> expanded_pars( max_num_par, 0.0 );
      vector<bool> is_used_parameter( max_num_par, false ), is_fit_parameter( max_num_par, false );
      size_t working_actual_index = 0; //indexes into `m_rel_eff_eqn_coefficients` and/or `m_rel_eff_eqn_covariance`
      
      // Fill out info on self-atten shielding
      if( is_valid_shield(m_phys_model_self_atten_shield.get(), m_input.phys_model_self_atten.get()) )
      {
        assert( !!m_input.phys_model_self_atten && !!m_phys_model_self_atten_shield );
        if( m_input.phys_model_self_atten->fit_atomic_number )
        {
          par_index[0] = working_actual_index;
          expanded_pars[0] = m_phys_model_self_atten_shield->m_atomic_number;
          is_used_parameter[0] = true;
          is_fit_parameter[0] = true;
          working_actual_index += 1;
        }
        
        is_used_parameter[1] = true;
        is_fit_parameter[1] = m_input.phys_model_self_atten->fit_areal_density;
        expanded_pars[1] = m_phys_model_self_atten_shield->m_areal_density;
        working_actual_index += 1;
      }//if( is_valid_shield(m_phys_model_self_atten_shield.get(), m_input.phys_model_self_atten.get()) )
      
      // Fill out info on external shieldings
      assert( m_input.phys_model_external_attens.size() == m_phys_model_external_atten_shields.size() );
      for( size_t i = 0; i < m_phys_model_external_atten_shields.size(); ++i )
      {
        const size_t expanded_an_index = 2 + 2*i;
        const size_t expanded_ad_index = expanded_an_index + 1;
        
        const auto &input = m_input.phys_model_external_attens[i];
        const auto &result = m_phys_model_external_atten_shields[i];
        assert( is_valid_shield(result.get(), input.get()) );
        if( !is_valid_shield(result.get(), input.get()) )
          continue;
        
        if( input->fit_atomic_number )
        {
          par_index[expanded_an_index] = working_actual_index;
          expanded_pars[expanded_an_index] = result->m_atomic_number;
          is_used_parameter[expanded_an_index] = true;
          is_fit_parameter[expanded_an_index] = true;
          working_actual_index += 1;
        }
        
        is_used_parameter[expanded_ad_index] = true;
        is_fit_parameter[expanded_ad_index] = input->fit_areal_density;
        expanded_pars[expanded_ad_index] = result->m_areal_density;
        working_actual_index += 1;
      }//for( size_t i = 0; i < m_phys_model_external_atten_shields.size(); ++i )
      
      // Fill out info on modified Hoerl function
      if( m_input.phys_model_use_hoerl )
      {
        assert( (working_actual_index + 1) < m_rel_eff_eqn_coefficients.size() );
        
        const size_t expanded_b_index = 2 + 2*m_phys_model_external_atten_shields.size();
        const size_t expanded_c_index = expanded_b_index + 1;
        assert( (expanded_c_index+1) == max_num_par );
        
        is_used_parameter[expanded_b_index] = true;
        is_fit_parameter[expanded_b_index] = true;
        expanded_pars[expanded_b_index] = m_rel_eff_eqn_coefficients.at(working_actual_index);
        par_index[expanded_b_index] = working_actual_index;
        working_actual_index += 1;
        
        
        is_used_parameter[expanded_c_index] = true;
        is_fit_parameter[expanded_c_index] = true;
        expanded_pars[expanded_c_index] = m_rel_eff_eqn_coefficients.at(working_actual_index);
        par_index[expanded_c_index] = working_actual_index;
        working_actual_index += 1;
      }//if( m_input.phys_model_use_hoerl )
      
      assert( working_actual_index == num_actual_pars );
      if( working_actual_index != num_actual_pars )
        throw std::logic_error( "eval_eqn_uncertainty:  working_actual_index != num_actual_pars" );
      
      const float energyf = static_cast<float>(energy);
      
      // Returns attenuation coefficient, mu, for expanded parameter number passed in
      //  lambda name corresponds to notation used to get the derivatives in the online derivative
      //  tool.
      const auto function_m = [&]( size_t expanded_index ) -> double {
        assert( (expanded_index % 2) == 0 );
        
        assert( (expanded_index == 0)
               || ((expanded_index - 2)/2) < m_phys_model_external_atten_shields.size() );
        if( (expanded_index != 0)
               && ((expanded_index - 2)/2) >= m_phys_model_external_atten_shields.size() )
          throw std::logic_error( "(expanded_index - 2)/2) >= m_phys_model_external_atten_shields.size()" );
          
        
        const auto &shield = (expanded_index == 0) ? m_phys_model_self_atten_shield
                                                     : m_phys_model_external_atten_shields[(expanded_index - 2)/2];
        assert( shield );
        if( !shield )
          throw std::logic_error( "logic error: !shield" );
        
        const auto &mat = shield->m_material;
        if( mat )
          return GammaInteractionCalc::transmition_length_coefficient( mat.get(), energyf ) / mat->density;
        return RelActCalc::get_atten_coef_for_an<double>( shield->m_atomic_number , energyf );
      };//const auto function_m
      
      const auto derivative_function_m = [&]( size_t expanded_index ) -> double {
        assert( (expanded_index == 0)
               || ((expanded_index - 2)/2) < m_phys_model_external_atten_shields.size() );
        if( (expanded_index != 0)
               && ((expanded_index - 2)/2) >= m_phys_model_external_atten_shields.size() )
          throw std::logic_error( "(expanded_index - 2)/2) >= m_phys_model_external_atten_shields.size()" );
        
        const auto &shield = (expanded_index == 0) ? m_phys_model_self_atten_shield
                                                     : m_phys_model_external_atten_shields[(expanded_index - 2)/2];
        assert( shield && !shield->m_material );
        
        ceres::Jet<double, 1> x(shield->m_atomic_number, 0); // x = AN, derivative index = 0
        ceres::Jet<double, 1> y = RelActCalc::get_atten_coef_for_an( x, energyf );
        //double value = y.a; // The function value at AN
        const double derivative = y.v[0]; // The derivative at AN
        
        //Printing things out, like below, confirms we are getting the correct derivative
        //cout << "Derivative for energy=" << energy << ", AN=" << shield->m_atomic_number << ", mu=" << y.a << " is: " << derivative << endl;
        //if( shield->m_atomic_number > 2 && shield->m_atomic_number < 96 )
        //{
        //  const int lower_an = std::max( 1, static_cast<int>( std::floor(shield->m_atomic_number) ) );
        //  const int upper_an = lower_an + 1;
        //  cout << "mu(" << lower_an << ")= " << MassAttenuation::massAttenuationCoeficient(lower_an, energy) 
        //       << ", mu(" << upper_an << ")= " << MassAttenuation::massAttenuationCoeficient(upper_an, energy) 
        //       << " (diff=" << (MassAttenuation::massAttenuationCoeficient(lower_an, energy) - MassAttenuation::massAttenuationCoeficient(upper_an, energy)) << ")" 
        //       << endl;
        //}
        
        return derivative;
      };//derivative_function_m(...)
      
      
      
      // Get derivative of f, with respect to parameter `expanded_index`.
      //  As of 20250114, totally untested
      auto df_dp = [&]( const size_t expanded_index ) -> double {
        assert( expanded_index < is_used_parameter.size() );
        assert( is_used_parameter[expanded_index] && is_fit_parameter[expanded_index] );
        
        const double det_eff = (m_input.phys_model_detector && m_input.phys_model_detector->isValid())
                                ? m_input.phys_model_detector->intrinsicEfficiency(energyf)
                                : 1.0f;
        
        double self_atten_val = 1.0;
        if( is_used_parameter[1] )
        {
          const double mu_0 = function_m(0);
          const double x_1 = expanded_pars[1];
          self_atten_val = ((1 - exp(-mu_0 * x_1)) / (mu_0 * x_1));
        }
        
        double hoerl_val = 1.0;
        if( m_input.phys_model_use_hoerl )
        {
          const size_t expanded_b_index = 2 + 2*m_phys_model_external_atten_shields.size();
          const size_t expanded_c_index = expanded_b_index + 1;
          const double b = expanded_pars[expanded_b_index];
          const double c = expanded_pars[expanded_c_index];
          hoerl_val = std::pow( (energy/1000.0), b) * std::pow( c, (1000.0 / energy) );
        }
        
        double ext_atten_exp_arg = 0.0;
        for( size_t ext_shield_ind = 0; ext_shield_ind < m_phys_model_external_atten_shields.size(); ++ext_shield_ind )
        {
          const size_t an_index = 2 + 2*ext_shield_ind;
          const size_t ad_index = an_index + 1;
          ext_atten_exp_arg -= function_m(an_index) * expanded_pars[ad_index];
        }
        const double ext_atten_val = exp( ext_atten_exp_arg );
        
        
        if( expanded_index == 0 )
        {
          // We want the derivative with respect to self-attenuating AN
          const double x_1 = expanded_pars[1]; //self-atten AD
          const double mu_0 = function_m(0);
          double exp_arg = -x_1 * mu_0;
          for( size_t ext_ind = 0; ext_ind < m_phys_model_external_atten_shields.size(); ++ext_ind )
          {
            const size_t expanded_an_index = 2 + 2*ext_ind;
            const size_t expanded_ad_index = expanded_an_index + 1;
            exp_arg -= function_m(expanded_an_index) * expanded_pars[expanded_ad_index];
          }
          
          double answer = det_eff * (exp(exp_arg) * (x_1 * mu_0 - exp(x_1 * mu_0) + 1) * derivative_function_m(0))*hoerl_val
                          / (x_1 * mu_0 * mu_0);
          
          return answer;
        }else if( expanded_index == 1 )
        {
          // We want the derivative with respect to self-attenuating AD
          const double x_1 = expanded_pars[1]; //self-atten AD
          const double mu_0 = function_m(0);
          
          double exp_arg = -mu_0 * x_1;
          for( size_t ext_ind = 0; ext_ind < m_phys_model_external_atten_shields.size(); ++ext_ind )
          {
            const size_t expanded_an_index = 2 + 2*ext_ind;
            const size_t expanded_ad_index = expanded_an_index + 1;
            exp_arg -= function_m(expanded_an_index) * expanded_pars[expanded_ad_index];
          }
      
          double answer = -det_eff * (exp(mu_0 * x_1) - mu_0 * x_1 - 1) * exp( exp_arg ) * hoerl_val
                          / (mu_0 * x_1*x_1);
          
          return answer;
        }else if( expanded_index < (2 + 2*m_phys_model_external_atten_shields.size()) )
        {
          const size_t ext_shield_num = (expanded_index - 2) / 2;
          
          double answer = -det_eff * hoerl_val * self_atten_val * ext_atten_val;
          
          if( ((expanded_index % 2) == 0) )
          {
            //External shield AN
            const double x_ad = expanded_pars[expanded_index + 1];
            answer *= x_ad * derivative_function_m(expanded_index);
          }else
          {
            //External shield AD
            answer *= function_m(expanded_index-1);
          }
          
          return answer;
        }
        
        //Modified Hoerl function
        if( ((expanded_index % 2) == 0) )
        {
          //df/dx_6 = (E^x_6 * (ln(E) - ln(1000)) * (e^(function_m(x_0) * x_1) - 1) * e^(-function_m(x_4) * x_5 - function_m(x_2) * x_3 - function_m(x_0) * x_1) * x_7^(1000 / E))
          //    / (function_m(x_0) * x_1 * 1000^x_6)
          
          double answer = det_eff * ext_atten_val * self_atten_val * hoerl_val * (log(energy) - log(1000));
          return answer;
        }
        
        const size_t expanded_b_index = 2 + 2*m_phys_model_external_atten_shields.size();
        const size_t expanded_c_index = expanded_b_index + 1;
        const double b = expanded_pars[expanded_b_index];
        const double c = expanded_pars[expanded_c_index];
          
        double answer = det_eff * ext_atten_val * self_atten_val * std::pow( (energy/1000.0), b);
        answer *= (1000.0/energy)*std::pow(c, (1000.0/energy) - 1.0 );
          
        return answer;
      };//auto df_dp = ...
      
      
      size_t i_working = 0;
      for( size_t i_expanded = 0; i_expanded < max_num_par; ++i_expanded )
      {
        if( !is_used_parameter[i_expanded] )
          continue;
        
        if( !is_fit_parameter[i_expanded] ) //We only include AN as a parameter when we fit it
        {
          assert( m_rel_eff_eqn_covariance[i_working][i_working] == 0.0 );
          
          i_working += 1;
          continue;
        }
        
        size_t j_working = 0;
        for( size_t j_expanded  = 0; j_expanded < max_num_par; ++j_expanded )
        {
          if( !is_used_parameter[j_expanded] )
            continue;
          
          if( !is_fit_parameter[j_expanded] ) //We only include AN as a parameter when we fit it
          {
            j_working += 1;
            continue;
          }
          
          const double i_component = df_dp( i_expanded );
          const double j_component = df_dp( j_expanded );
          
          uncert_sq += i_component * m_rel_eff_eqn_covariance[i_working][j_working] * j_component;
          
          j_working += 1;
        }//for( size_t j_expanded  = 0; j_expanded < max_num_par; ++j_expanded )
        assert( j_working == num_actual_pars );
        
        i_working += 1;
      }//for( size_t i_expanded  = 0; i_expanded < max_num_par; ++i_expanded )
      assert( i_working == num_actual_pars );
      
      break;
    }//case RelActCalc::RelEffEqnForm::FramPhysicalModel:
  }//switch( eqn_form )

  const double eval_val = rel_eff_eqn_value( energy );
  //cout << "Phys Model Rel Eff: energy=" << energy << "keV, val=" << eval_val << " uncert_sq=" << uncert_sq << " --> uncert=" << sqrt(uncert_sq) << endl;

  //assert( uncert_sq >= 0.0 );
  if( uncert_sq < 0.0 )
    throw std::runtime_error( "RelEffSolution::rel_eff_eqn_uncert: negative squared uncertainty." );
  
  return sqrt(uncert_sq);
}//double rel_eff_eqn_uncert( const double energy ) const


string RelEffSolution::rel_eff_eqn_js_uncert_fcn() const
{
  vector<double> energies;
  double current_energy = 20.0;
  const double upper_energy = m_input.peaks.empty() ? 3000.0 
                                                    : std::max(3000.0, m_input.peaks.back().m_energy);

  for( const GenericPeakInfo &peak : m_input.peaks )
  {
    double min_dx = 1.0;
    if( current_energy < 130 )
      min_dx = 1.0;
    else if( current_energy < 300 )
      min_dx = 5.0;
    else
      min_dx = 15.0;

    // We'll try to get in at least ~10 points between each peak
    if( peak.m_energy > current_energy )
      min_dx = std::min( min_dx, 0.1*(peak.m_energy - current_energy) );
    min_dx = std::max( min_dx, 1.0 ); //but less than a keV between points is just to small.

    for( ; current_energy < peak.m_energy; current_energy += min_dx )
      energies.push_back( current_energy );
    
    if( !energies.empty() && (energies.back() < peak.m_energy) )
      current_energy = peak.m_energy;
  }//for( const GenericPeakInfo &peak : m_input.peaks )
  
  for( ; current_energy < upper_energy; current_energy += 15 )
    energies.push_back( current_energy );
  
  size_t num_points = 0;
  string fcn = "function(x){\n"
  "  const points = [";
  bool is_first_point = true;
  for( double x : energies )
  {
    try
    {
      double y = rel_eff_eqn_uncert( x );
      //assert( (y >= 0.0) && !IsNan(y) && !IsInf(y) ); //Can happen when we are out of bounds of the physical model
      if( isnan(y) || isinf(y) )
        continue;
    
      fcn += is_first_point ? "" : ",";
      fcn += "[" + SpecUtils::printCompact(x, 4) + "," + SpecUtils::printCompact(y, 4) + "]";

      num_points += 1;
      is_first_point = false;
    }catch( std::exception &e )
    {
      // This can happen when we are out of bounds of the physical model
    }
  }//for( double x : energies )

if( num_points < 2 )
 return "null";

  fcn += "];\n"
  "  if( x <= points[0][0] )\n"
  "    return points[0][1];\n"
  "  if( x >= points[points.length - 1][0] )\n"
  "    return points[points.length - 1][1];\n"
  "  for (let i = 0; i < points.length - 1; i++) {\n"
  "    const [x1, y1] = points[i];\n"
  "    const [x2, y2] = points[i + 1];\n"
  "    if( x >= x1 && x <= x2) {\n"
  "      const t = (x - x1) / (x2 - x1);\n"
  "      return y1 + t * (y2 - y1);\n"
  "    }\n"
  "  }\n"
  "console.assert(0,'Shouldnt get here in interpolating');\n"
  "return points[points.length - 1][1];"
  "}";
  
  return fcn;
}//string RelEffSolution::rel_eff_eqn_js_uncert_fcn() const


string RelEffSolution::rel_eff_eqn_txt( const bool html_format ) const
{
  if( m_input.eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    const ManualGenericRelActFunctor::PhysModelRelEqnDef input
          = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );
    return RelActCalc::physical_model_rel_eff_eqn_text( input.self_atten,
                input.external_attens, input.det, input.hoerl_b, input.hoerl_b, html_format );
  }
  
  return RelActCalc::rel_eff_eqn_text( m_input.eqn_form, m_rel_eff_eqn_coefficients );
}//std::string RelEffSolution::rel_eff_eqn_txt( const bool html_format ) const
  
  
string RelEffSolution::rel_eff_eqn_js_function() const
{
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
    return RelActCalc::rel_eff_eqn_js_function( m_input.eqn_form, m_rel_eff_eqn_coefficients );
  

  const ManualGenericRelActFunctor::PhysModelRelEqnDef input
        = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );
    
  return RelActCalc::physical_model_rel_eff_eqn_js_function( input.self_atten,
                          input.external_attens, input.det.get(), input.hoerl_b, input.hoerl_c );
}//string RelEffSolution::rel_eff_eqn_js_function() const
  

void RelEffSolution::print_html_report( ostream &output_html_file,
                                       string spectrum_title,
                                       shared_ptr<const SpecUtils::Measurement> spectrum,
                                       vector<shared_ptr<const PeakDef>> display_peaks,
                                       shared_ptr<const SpecUtils::Measurement> background,
                                       double background_normalization ) const
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  const size_t nsigfig = 4;
  
  char buffer[512] = { '\0' };
  
  
  stringstream results_html;
  results_html << "<div>&chi;<sup>2</sup>=" << SpecUtils::printCompact(m_chi2, nsigfig)
  << " and there were " << m_dof << " DOF; &chi;<sup>2</sup>/dof="
  << SpecUtils::printCompact(m_chi2/m_dof, nsigfig)
  << " </div>\n";
  
  results_html << "<div class=\"releffeqn\">Rel. Eff. Eqn: y = ";
  results_html << rel_eff_eqn_txt( true );
  
  results_html << "</div>\n";
  
  get_mass_fraction_table( results_html );
  get_mass_ratio_table( results_html );
  
  
  const bool has_decay_corr = !m_input.peaks_before_decay_correction.empty();
  
  // Make table giving info on each of the _used_ peaks
  results_html << "<table class=\"peaktable resulttable\">\n";
  results_html << "  <caption>Peaks used for analysis.</caption>\n";
  results_html << "  <thead><tr>"
  "<th scope=\"col\">Energy (keV)</th>"
  "<th scope=\"col\">Nuclide</th>"
  "<th scope=\"col\">Yield</th>"
  "<th scope=\"col\">Net Area</th>"
  "<th scope=\"col\">Net Area Uncert</th>"
  "<th scope=\"col\">Counts/Yield</th>"
  "<th scope=\"col\">Counts/Yield Unc</th>"
  "<th scope=\"col\">Add. Unc.</th>"
  "<th scope=\"col\">Meas. Rel Eff</th>"
  "<th scope=\"col\">Meas. Rel Eff Unct</th>"
  << (has_decay_corr ? "<th scope=\"col\">Decay Corr.</th>" : "")
  << "</tr></thead>\n"
  "  <tbody>\n";
  
  
  for( const GenericPeakInfo &info : m_input.peaks )
  {
    snprintf(buffer, sizeof(buffer), "%.2f", info.m_mean );
    results_html << "  <tr><td>" << buffer << "</td>";
    for( size_t i = 0; i < info.m_source_gammas.size(); ++i )
    {
      const GenericLineInfo &line = info.m_source_gammas[i];
      
      if( i )
        results_html << "<tr><td></td>";
      
      const double rel_act = relative_activity(line.m_isotope);
      const double counts_over_yield = info.m_counts / line.m_yield;
      const double counts_uncert_percent = 100.0*(info.m_counts_uncert / info.m_counts);
      const double meas_rel_eff = info.m_counts / (line.m_yield * rel_act);
      const double meas_rel_eff_uncert = 100* info.m_counts_uncert / info.m_counts;
      
      results_html << "<td>" << line.m_isotope
      << "</td><td>" << SpecUtils::printCompact( line.m_yield, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( info.m_counts, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( info.m_counts_uncert, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( counts_over_yield, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( counts_uncert_percent, nsigfig ) << "%"
      << "</td><td>" << SpecUtils::printCompact( info.m_base_rel_eff_uncert, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( meas_rel_eff, nsigfig )
      << "</td><td>" << SpecUtils::printCompact( meas_rel_eff_uncert, nsigfig ) << "%";
      
      if( has_decay_corr )
      {
        results_html << "</td><td>";
        
        for( const GenericPeakInfo &un_corr_peak : m_input.peaks_before_decay_correction )
        {
          if( fabs(info.m_energy - un_corr_peak.m_energy) > 0.001 )
            continue;
          
          for( const GenericLineInfo &un_corr_line : un_corr_peak.m_source_gammas )
          {
            if( un_corr_line.m_isotope != line.m_isotope )
              continue;
            
            const double ratio = line.m_yield / un_corr_line.m_yield;
            if( !IsInf(ratio) && !IsNan(ratio) )
              results_html << " " << SpecUtils::printCompact( ratio, nsigfig );
            else
              results_html << "--";
          }//for( loop over un_corr_peak.m_source_gammas )
        }//for( loop over m_input.peaks_before_decay_correction )
      }//if( has_decay_corr )
      
      results_html << "</td></tr>\n";
    }//for( size_t i = 0; i < info.m_source_gammas.size(); ++i )
  }//for( const RelEff::PeakInfo &info : input_peaks )
  
  results_html << "  </tbody>\n"
  << "</table>\n\n";
  
  
  auto html_sanitize = []( string &val ){
    // We'll do some really stupid/simple sanitization
    SpecUtils::ireplace_all( val, "&", "&amp;" );
    SpecUtils::ireplace_all( val, "<", "&lt;" );
    SpecUtils::ireplace_all( val, ">", "&gt;" );
    SpecUtils::ireplace_all( val, "'", "&#39;" );
    SpecUtils::ireplace_all( val, "\"", "&quot;" );
  };
  
  results_html << "</div>\n";
  
  
  if( !m_warnings.empty() )
  {
    results_html << "<div class=\"warnings\">\n"
    << "<h3>Warnings</h3>\n";
    for( string warning : m_warnings )
    {
      html_sanitize( warning );
      
      results_html << "<div class=\"warningline\">" << warning << "</div>\n";
    }//for( string warning : warnings )
    
    results_html << "</div>\n";
  }//if( !warnings.empty() )
  
  
  time_t rawtime;
  struct tm *timeinfo;
  time( &rawtime );
  timeinfo = localtime (&rawtime);
  results_html << "<div class=\"anatime\">Analysis performed " << asctime(timeinfo)
  << " with " << spectrum_title << " compiled " __TIMESTAMP__
  << "</div>\n";
  
  
  //Write out the data JSON
  stringstream rel_eff_plot_values, add_rel_eff_plot_css;
  rel_eff_plot_values << "[";
  for( size_t index = 0; index < m_input.peaks.size(); ++index )
  {
    const GenericPeakInfo &peak = m_input.peaks[index];
    
    string isotopes_json;
    double src_counts = 0.0;
    for( const GenericLineInfo &line : peak.m_source_gammas )
    {
      const double rel_act = relative_activity( line.m_isotope );
      src_counts += rel_act * line.m_yield;
      
      //const double meas_rel_eff = info.m_counts / (info.m_source_gammas[i].m_yield * rel_act);
      
      snprintf( buffer, sizeof(buffer), "%s{\"nuc\": \"%s\", \"br\": %1.6G, \"rel_act\": %1.6G}",
               (isotopes_json.empty() ? "" : ", "), line.m_isotope.c_str(), line.m_yield, rel_act );
      
      isotopes_json += buffer;
    }//for( const RelEff::GammaLineInfo &line : peak.m_source_gammas )
    
    const double eff = peak.m_counts / src_counts;
    const double eff_uncert = peak.m_counts_uncert / src_counts;
    
    snprintf( buffer, sizeof(buffer),
             "%s"
             "{\"energy\": %.2f, \"counts\": %1.7g,"
             " \"counts_uncert\": %1.7g, \"eff\": %1.6g,"
             " \"eff_uncert\": %1.6g",
             (index ? ", " : ""), 
             peak.m_mean, peak.m_counts,
             peak.m_counts_uncert, eff,
             eff_uncert );
    
    rel_eff_plot_values << buffer;
    rel_eff_plot_values << ", \"nuc_info\": [" << isotopes_json << "]}";
  }//for( size_t index = 0; index < input_peaks.size(); ++index )
  
  rel_eff_plot_values << "]";
  
  
  set<const SandiaDecay::Nuclide *> nuclides_with_colors;
  for( const shared_ptr<const PeakDef> &p : display_peaks )
  {
    assert( p );
    const SandiaDecay::Nuclide * const nuc = p ? p->parentNuclide() : nullptr;
    if( nuc && !nuclides_with_colors.count(nuc) && !p->lineColor().isDefault() )
    {
      add_rel_eff_plot_css << "        .RelEffPlot circle." << nuc->symbol
                           << "{ fill: " << p->lineColor().cssText() << "; }\n";
      nuclides_with_colors.insert( nuc );
    }
  }//for( const shared_ptr<const PeakDef> &p : display_peaks )
  
  size_t unseen_nuc_index = 0;
  const vector<string> default_nuc_colors{ "#003f5c", "#ffa600", "#7a5195", "#ff764a", "#ef5675", "#374c80" };
  for( const IsotopeRelativeActivity &act : m_rel_activities )
  {
    const SandiaDecay::Nuclide * const nuc = db->nuclide( act.m_isotope );
    if( !nuclides_with_colors.count(nuc) )
    {
      string nucstr = act.m_isotope;
      if( nucstr.empty() )
        nucstr = "default";
      
      for( size_t i = 0; i < nucstr.size(); ++i )
        nucstr[i] = std::isalpha( static_cast<unsigned char>(nucstr[i]) ) ? nucstr[i] : '_';
      
      add_rel_eff_plot_css << "        .RelEffPlot circle." << nucstr << "{ fill: "
      << default_nuc_colors[unseen_nuc_index % default_nuc_colors.size()]
      << "; }\n";
      
      unseen_nuc_index += 1;
    }
  }//for( const IsotopeRelativeActivity &act : m_rel_activities )
  
  
  auto load_file_contents = []( string filename ) -> string {
    Wt::WApplication *app = Wt::WApplication::instance();
    string filepath = app ? app->docRoot() : string("");
    filepath = SpecUtils::append_path( filepath, "InterSpec_resources" );
    filepath = SpecUtils::append_path( filepath, filename );
    
    vector<char> file_data;
    try
    {
      SpecUtils::load_file_data( filepath.c_str(), file_data );
    }catch( std::exception & )
    {
      throw std::runtime_error( "Failed to read " + filename );
    }
    
    return string( begin( file_data ), end( file_data ) );
  };//load_file_contents(...)
  
  
  string html = load_file_contents( "static_text/manual_rel_act_report.tmplt.html" );
  const string d3_js = load_file_contents( "d3.v3.min.js" );
  
  SpecUtils::ireplace_all( html, "\\;", ";" );
  
  SpecUtils::ireplace_all( html, "${D3_SCRIPT}", d3_js.c_str() );
  
  SpecUtils::ireplace_all( html, "${TITLE}", spectrum_title.c_str() );
  
  SpecUtils::ireplace_all( html, "${REL_EFF_DATA_VALS}", rel_eff_plot_values.str().c_str() );
  
  string rel_eff_fcn;
  if( m_input.eqn_form != RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    rel_eff_fcn = RelActCalc::rel_eff_eqn_js_function( m_input.eqn_form, m_rel_eff_eqn_coefficients );
  }else
  {
    const ManualGenericRelActFunctor::PhysModelRelEqnDef input
           = ManualGenericRelActFunctor::make_phys_eqn_input( m_input, m_rel_eff_eqn_coefficients );
    
    rel_eff_fcn = RelActCalc::physical_model_rel_eff_eqn_js_function( input.self_atten,
                                                                     input.external_attens,
                                                                     input.det.get(),
                                                                     input.hoerl_b,
                                                                     input.hoerl_c );
  }

  string unc_fcn = rel_eff_eqn_js_uncert_fcn();

  SpecUtils::ireplace_all( html, "${FIT_REL_EFF_EQUATION}", rel_eff_fcn.c_str() );
  SpecUtils::ireplace_all( html, "${FIT_REL_EFF_EQUATION_UNCERTAINTY}", unc_fcn.c_str() );
  SpecUtils::ireplace_all( html, "${RESULTS_TXT}", results_html.str().c_str() );
  
  
  const string rel_eff_plot_js = load_file_contents( "RelEffPlot.js" );
  const string rel_eff_plot_css = load_file_contents( "RelEffPlot.css" );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_JS}", rel_eff_plot_js.c_str() );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_CSS}", rel_eff_plot_css.c_str() );
  SpecUtils::ireplace_all( html, "${REL_EFF_PLOT_ADDITIONAL_CSS}", add_rel_eff_plot_css.str().c_str() );
  
  
  if( spectrum )
  {
    stringstream set_js_str;
    
    D3SpectrumExport::write_js_for_chart( set_js_str, "specchart", "", "Energy (keV)", "Counts"  );
    
    set_js_str <<
    "  let spec_observer = new ResizeObserver(entries => {\n"
    "    for (let entry of entries) {\n"
    "      if (entry.target && (entry.target.id === \"specchart\")) {\n"
    "        spec_chart_specchart.handleResize(false);\n"
    "      }\n"
    "    }\n"
    "  });\n"
    "  spec_observer.observe( document.getElementById(\"specchart\") );\n"
    ;
    
    D3SpectrumExport::D3SpectrumChartOptions chart_options;
    chart_options.m_useLogYAxis = true;
    chart_options.m_legendEnabled = false;
    chart_options.m_compactXAxis = true;
    chart_options.m_allowDragRoiExtent = false;
    write_set_options_for_chart( set_js_str, "specchart", chart_options );
    set_js_str << "  spec_chart_specchart.setShowLegend(false);\n";
    set_js_str << "  spec_chart_specchart.setXAxisRange("
               << spectrum->gamma_energy_min()
               << ", " << spectrum->gamma_energy_max() << ", false);\n";
    
    D3SpectrumExport::D3SpectrumOptions spec_options;
    spec_options.spectrum_type = SpecUtils::SpectrumType::Foreground;
    
    vector<shared_ptr<const PeakDef>> peaks;
    for( const auto &p : display_peaks )
      peaks.push_back( make_shared<PeakDef>(*p) );
    spec_options.peaks_json = PeakDef::peak_json( peaks, spectrum );
    
    
    vector<pair<const SpecUtils::Measurement *,D3SpectrumExport::D3SpectrumOptions> > meas_to_plot;
    
    const SpecUtils::Measurement * const meas_ptr = spectrum.get();
    meas_to_plot.emplace_back(meas_ptr, spec_options);
    
    if( background )
    {
      D3SpectrumExport::D3SpectrumOptions spec_options;
      spec_options.spectrum_type = SpecUtils::SpectrumType::Background;
      spec_options.line_color = "steelblue";
      spec_options.display_scale_factor = 1.0;
      
      if( (background_normalization <= 0.0f)
         || IsNan(background_normalization)
         || IsInf(background_normalization) )
      {
        float back_lt = background->live_time();
        if( (back_lt <= 0.0f) || IsNan(back_lt) || IsNan(back_lt) )
          back_lt = background->real_time();
        
        float fore_lt = spectrum->live_time();
        if( (fore_lt <= 0.0f) || IsNan(fore_lt) || IsNan(fore_lt) )
          fore_lt = spectrum->real_time();
        
        background_normalization = fore_lt / back_lt;
      }//if( background normalization wasnt provided )
      
      if( (background_normalization > 0.0f)
         && !IsNan(background_normalization)
         && !IsInf(background_normalization) )
      {
        spec_options.display_scale_factor = background_normalization;
      }
        
      spec_options.title = "Background";
      meas_to_plot.emplace_back(background.get(), spec_options);
    }//if( background )
    
    
    D3SpectrumExport::write_and_set_data_for_chart( set_js_str, "specchart", meas_to_plot );
    
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_DIV}",
                            "<div id=\"specchart\" style=\"height: 30vw; flex: 1 2; overflow: hidden;\" class=\"SpecChart\"></div>" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_INIT_JS}", set_js_str.str().c_str() );
    
    
    const string spectrum_chart_d3_js = load_file_contents( "SpectrumChartD3.js" );
    const string spectrum_chart_d3_css = load_file_contents( "SpectrumChartD3.css" );
    
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_JS}", spectrum_chart_d3_js.c_str() );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_CSS}", spectrum_chart_d3_css.c_str() );
    
    SpecUtils::ireplace_all( html, "${CHART_SPACER_LEFT}", "" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_RIGHT}", "" );
  }else
  {
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_DIV}", "" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_JS}", "" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_CSS}", "" );
    SpecUtils::ireplace_all( html, "${SPECTRUM_CHART_INIT_JS}", "" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_LEFT}", "<div style=\"width: 10%\"> </div>" );
    SpecUtils::ireplace_all( html, "${CHART_SPACER_RIGHT}", "<div style=\"width: 15%\"> </div>" );
  }//if( spectrum ) / else
  
  
  output_html_file << html;
}//void print_html_report( std::ostream &strm ) const


void setup_physical_model_shield_par( ceres::Problem * const problem, 
                                      double * const pars, 
                                      const size_t start_ind, 
                                      const std::shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt )
{
  const size_t an_index = start_ind + 0;
  const size_t ad_index = start_ind + 1;
  
  if( !opt || (!opt->material && (opt->atomic_number == 0.0)) )
  {
    pars[an_index] = 0.0;
    pars[ad_index] = 0.0;
    if( problem )
    {
      problem->SetParameterBlockConstant( pars + an_index );
      problem->SetParameterBlockConstant( pars + ad_index );
    }
    return;
  }//if( !opt )

  opt->check_valid();

  if( opt->material )
  {
    if( opt->fit_atomic_number )
      throw runtime_error( "You can not fit AN when defining a material" );

    pars[an_index] = 0.0;       
    if( problem )
      problem->SetParameterBlockConstant( pars + an_index );
  }else
  {     
    double lower_an = opt->lower_fit_atomic_number;
    double upper_an = opt->upper_fit_atomic_number;
    if( (lower_an == upper_an) && (lower_an == 0.0) )
    {
      lower_an = 1.0;
      upper_an = 98.0;
    }

    double an = opt->atomic_number;
    if( (an == 0.0) && opt->fit_atomic_number )
      an = 0.5*(opt->lower_fit_atomic_number + opt->upper_fit_atomic_number);
          
    if( (an < 1) || (an > 98) || (an < lower_an) || (an > upper_an) )
      throw runtime_error( "Self-atten AN is invalid" );
          
    pars[an_index] = an / RelActCalc::ns_an_ceres_mult;
          
    if( opt->fit_atomic_number )
    {
      if( (lower_an < 1) || (upper_an > 98) || (upper_an <= lower_an) )
        throw runtime_error( "Self-atten AN limits is invalid" );
            
      if( problem )
      {
        problem->SetParameterLowerBound( pars + an_index, 0, lower_an / RelActCalc::ns_an_ceres_mult );
        problem->SetParameterUpperBound( pars + an_index, 0, upper_an / RelActCalc::ns_an_ceres_mult );

        // We could Quantize AN... but I'm not sure the below does it
        //std::vector<int> quantized_params(1, 0);  // 0 means continuous, 1 means quantized
        //quantized_params[an_index] = 1;  // Assuming an_index is the index of the AN parameter
        //problem->SetParameterization( pars + an_index, 
        //     new ceres::SubsetParameterization(1, quantized_params));
      }
    }else
    {
      if( problem )
        problem->SetParameterBlockConstant( pars + an_index );
    }
  }//if( opt.material ) / else
      
  const double max_ad = RelActCalc::PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2;
  double ad = opt->areal_density / PhysicalUnits::g_per_cm2;
  double lower_ad = opt->lower_fit_areal_density / PhysicalUnits::g_per_cm2;
  double upper_ad = opt->upper_fit_areal_density / PhysicalUnits::g_per_cm2;
        
  if( (lower_ad == upper_ad) && (lower_ad == 0.0) )
  {
    lower_ad = 0.0;
    upper_ad = max_ad;
  }
        
  if( (ad == 0.0) && opt->fit_areal_density )
  {
    ad = 2.5; // We want something away from zero, because Ceres doesnt like zero values much - 2.5 is arbitrary
    //  ad = 0.5*(lower_ad + upper_ad); //Something like 250 would be way too much
  }
  if( (ad < 0.0) || (ad > max_ad) )
    throw runtime_error( "Self-atten AD is invalid" );
        
  pars[ad_index] = ad;
        
  if( opt->fit_areal_density )
  {
    // Check for limits
    if( (lower_ad < 0.0) || (upper_ad > max_ad) || (lower_ad >= upper_ad) )
      throw runtime_error( "Self-atten AD limits is invalid" );
          
    if( problem )
    {
      problem->SetParameterLowerBound( pars + ad_index, 0, lower_ad );
      problem->SetParameterUpperBound( pars + ad_index, 0, upper_ad );
    }
  }else
  {
    if( problem )
      problem->SetParameterBlockConstant( pars + ad_index );
  }
}//void setup_physical_model_shield_par( ceres::Problem... )



RelEffSolution solve_relative_efficiency( const RelEffInput &input_orig )
{
  // When fitting the AN using the Physical Model, we can easily get caught in a local-minimum, and also
  // we dont currently have great control over the step sizes for AN/AD, so here is 
  // a work-around, to get a decent starting point for AN.
  // Right now we are scanning on AN, but scanning on AD as well would be better.
  // I assume that we could avoid this with a proper implementation of DynamicCostFunction.
#define SCAN_AN_FOR_BEST_FIT 1


#if( SCAN_AN_FOR_BEST_FIT )
  RelEffInput input = input_orig; //tmp copy for trying hack
#else
  const RelEffInput &input = input_orig;
#endif

  const std::vector<GenericPeakInfo> &peak_infos = input.peaks;
  const RelActCalc::RelEffEqnForm eqn_form = input.eqn_form;
  const size_t eqn_order = input.eqn_order;

  const auto start_time = std::chrono::high_resolution_clock::now();
  
  RelEffSolution solution;
  
  DoWorkOnDestruct setFinalTime( [&solution,start_time](){
    const auto end_time = std::chrono::high_resolution_clock::now();
    solution.m_num_microseconds_eval = std::chrono::duration<double, std::micro>(end_time - start_time).count();
  });
  
  solution.m_input = input;
  solution.m_status = ManualSolutionStatus::NotInitialized;

  solution.m_warnings.insert( begin(solution.m_warnings),
                      begin(input.prep_warnings), end(input.prep_warnings) );

#if( SCAN_AN_FOR_BEST_FIT )
  const bool scan_an_for_best_fit = (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel 
                                      && input.phys_model_self_atten
                                      && input.phys_model_self_atten->fit_atomic_number );
  if( scan_an_for_best_fit )
  {
    SpecUtilsAsync::ThreadPool threadpool;
    std::mutex best_chi2_mutex;
    double best_chi2 = std::numeric_limits<double>::max();
    double best_an = 0.0, best_ad = 0.0;
    
    const auto do_work = [&input, &best_chi2, &best_an, &best_ad, &best_chi2_mutex]( const double an ){
      RelEffInput new_input = input;
      auto new_self_atten = std::make_shared<RelActCalc::PhysicalModelShieldInput>(*input.phys_model_self_atten);
      new_self_atten->fit_atomic_number = false;
      new_self_atten->atomic_number = an;
      new_input.phys_model_self_atten = new_self_atten;
      RelEffSolution solution = solve_relative_efficiency( new_input );
      
      std::lock_guard<std::mutex> lock(best_chi2_mutex);
      
      if( (solution.m_status==ManualSolutionStatus::Success) && (solution.m_chi2 < best_chi2) )
      {
        best_chi2 = solution.m_chi2;
        best_an = an;
        assert( solution.m_phys_model_self_atten_shield );
        best_ad = solution.m_phys_model_self_atten_shield->m_areal_density / PhysicalUnits::g_per_cm2;
      }
    };//do_work

    double an_step = 5.0;
    double min_an = input.phys_model_self_atten->lower_fit_atomic_number;
    double max_an = input.phys_model_self_atten->upper_fit_atomic_number;
    
    for( double an = min_an; an <= max_an; an += an_step )
      threadpool.post( [&do_work, an](){ do_work(an); } );
    threadpool.join();

    if( best_chi2 != std::numeric_limits<double>::max() )
    {
      cout << "Initial best AN = " << best_an << " AD = " << best_ad << endl;
      best_chi2 = std::numeric_limits<double>::max();
      min_an = std::max( min_an, best_an - an_step );
      max_an = std::min( max_an, best_an + an_step );
      an_step = 1.0;
      for( double an = min_an; an <= max_an; an += an_step )
        threadpool.post( [&do_work, an](){ do_work(an); } );
      threadpool.join();
    }
    
    if( best_chi2 != std::numeric_limits<double>::max() )
    {
      cout << "Final best AN = " << best_an << " AD = " << best_ad << endl;
      auto new_self_atten = std::make_shared<RelActCalc::PhysicalModelShieldInput>(*input.phys_model_self_atten);
      //new_self_atten->fit_atomic_number = input.use_ceres_to_fit_eqn ? false : true;
      new_self_atten->atomic_number = best_an;
      new_self_atten->areal_density = best_ad * PhysicalUnits::g_per_cm2;
      input.phys_model_self_atten = new_self_atten;
      
      // // How many parameters, and how they are defined will change based on if we are fitting
      // //  atomic number - so we need to update the input stored in the solution to reflect what
      // //  we actually used.
      solution.m_input = input;
    }//if( best_chi2 != std::numeric_limits<double>::max() )
  }//if( scan_an_for_best_fit )
#endif //SCAN_AN_FOR_BEST_FIT

  try
  {
    if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {  
      if( !input.phys_model_detector || !input.phys_model_detector->isValid() )
        throw runtime_error( "You must specify a detector for the Physical Model." );
      
      // TODO: instead look at number of nuclides, and number of fit parameters, and use that to decide
      //       if we are okay - because maybe we want to use a single peak to get a curve, but not
      //       actually fit anything
      if( peak_infos.size() < 2 )
        throw runtime_error( "You must specify at least two peaks for the Physical Model." );

      if( !input.use_ceres_to_fit_eqn )
        throw logic_error( "You must specify to use Ceres to fit Rel. Eff. equation for the Physical Model." );
      
      //if( input.use_ceres_to_fit_eqn && input.phys_model_self_atten && input.phys_model_self_atten->fit_atomic_number )
      //  throw logic_error( "You can not use Ceres to fit atomic number of Physical relative efficiency equation." );
      
      if( input.phys_model_self_atten && !input.phys_model_self_atten->material
         && ((input.phys_model_self_atten->atomic_number < 0.999) || (input.phys_model_self_atten->atomic_number > 98.001)) )
        throw logic_error( "A self-attenuator is specified, but with no material or atomic number." );
      
      // TODO: add checks for external attenuators
    }else
    {
      if( input.phys_model_self_atten || !input.phys_model_external_attens.empty() )
        throw runtime_error( "Attenuations can only be specified for FramPhysicalModel." );
    }//if( RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorInitializing;
    solution.m_error_message = e.what();
    return solution;
  }//try / catch to check input


  ManualGenericRelActFunctor *cost_functor = nullptr;
  try
  {
    cost_functor = new ManualGenericRelActFunctor( input );
    solution.m_status = ManualSolutionStatus::ErrorFindingSolution;
    
    solution.m_warnings.insert( end(solution.m_warnings),
                               begin(cost_functor->m_setup_warnings),
                               end(cost_functor->m_setup_warnings) );
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorInitializing;
    solution.m_error_message = e.what();
    
    return solution;
  }//try / catch to setup the problem
  
  const size_t num_peaks = cost_functor->m_input.peaks.size();
  const size_t num_nuclides = cost_functor->m_isotopes.size();
  const size_t num_parameters = cost_functor->num_parameters();

  solution.m_activity_norms = cost_functor->m_rel_act_norms;

  // Relative activities multiples start out as 1.0 because ManualGenericRelActFunctor constructor
  //   estimates the activities for a flat rel eff = 1.0; see
  //   #ManualGenericRelActFunctor::m_rel_act_norms.
  vector<double> parameters( num_parameters, 1.0 );
  double *pars = &parameters[0];


  ceres::CostFunction *cost_function = nullptr;

  // From a few example cases inspected by hand, it looks like auto and numerical differentiation
  //  get the same answers, but auto diff is faster, and requires a lot fewer evaluations, so
  //  we'll just always use auto differentiation.
  const bool use_auto_diff = true;
  if( use_auto_diff )
  {
    ceres::DynamicAutoDiffCostFunction<ManualGenericRelActFunctor> *dyn_auto_diff_cost_function
          = new ceres::DynamicAutoDiffCostFunction<ManualGenericRelActFunctor>( cost_functor,
                                                                          ceres::TAKE_OWNERSHIP );
    // The number of residuals is the number of peaks, unless USE_RESIDUAL_TO_BREAK_DEGENERACY then
    //  we add one more residual to clamp the relative efficiency curve to 1.0 at the lowest energy.
 
    dyn_auto_diff_cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );
    dyn_auto_diff_cost_function->AddParameterBlock( static_cast<int>(num_parameters) );

    cost_function = dyn_auto_diff_cost_function;
  }else
  {
    ceres::NumericDiffOptions num_diff_options;

#if( SCAN_AN_FOR_BEST_FIT )
    if( scan_an_for_best_fit )
    {
      // TODO: need to more closely evaluate the step sizes
      num_diff_options.ridders_relative_initial_step_size = 0.1;
      num_diff_options.relative_step_size = 1.0E-3;
    }
#endif //SCAN_AN_FOR_BEST_FIT

    ceres::DynamicNumericDiffCostFunction<ManualGenericRelActFunctor> *dyn_num_diff_cost_function
          = new ceres::DynamicNumericDiffCostFunction<ManualGenericRelActFunctor>( cost_functor,
                                                        ceres::TAKE_OWNERSHIP, num_diff_options );
    dyn_num_diff_cost_function->SetNumResiduals( static_cast<int>(cost_functor->number_residuals()) );        
    dyn_num_diff_cost_function->AddParameterBlock( static_cast<int>(num_parameters) );     

    cost_function = dyn_num_diff_cost_function;
  }//if( use_auto_diff ) / else
    
  
  ceres::Problem problem;
  
  // TODO: investigate using a LossFunction - probably really need it
  //       The Huber and Cauchy functions dont seem to help a  ton on a single problem; probably need to define our own, probably based on Huber
  ceres::LossFunction *lossfcn = nullptr;
  //ceres::LossFunction *lossfcn = new ceres::HuberLoss(3.5);
  //ceres::LossFunction *lossfcn = new ceres::CauchyLoss(1.0);

  vector<double *> parameter_blocks;
  parameter_blocks.push_back( pars );
  problem.AddResidualBlock( cost_function, lossfcn, parameter_blocks[0] );
  
  problem.AddParameterBlock( pars, static_cast<int>(num_parameters) );
  // We dont seem to need to add paramater block to the problem - it picks it up from the cost function
  //problem.AddParameterBlock( pars, static_cast<int>(num_nuclides) );
  //for( size_t i = 0; i < num_phys_model_rel_eff_params; ++i )
  //  problem.AddParameterBlock( pars + num_nuclides + i, 1 );


  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    assert( input.use_ceres_to_fit_eqn );
    
    size_t par_num = num_nuclides;
    assert( par_num <= num_parameters );
    
    vector<int> constant_parameters;
    
    const auto &self_atten_opt = input.phys_model_self_atten;
    setup_physical_model_shield_par_manual( constant_parameters, pars, par_num, self_atten_opt );
    
    assert( par_num <= num_parameters );
    
    for( size_t ext_ind = 0; ext_ind < input.phys_model_external_attens.size(); ++ext_ind )
    {
      assert( par_num <= num_parameters );
      const auto &opt = input.phys_model_external_attens[ext_ind];
      setup_physical_model_shield_par_manual( constant_parameters, pars, par_num, opt );
    }//for( size_t ext_ind = 0; ext_ind < options.phys_model_external_atten.size(); ++ext_ind )

    assert( par_num <= num_parameters );
    
    if( input.phys_model_use_hoerl )
    {
      // set the b and c parameters for the relative efficiency equation
      pars[par_num] = 0.0;  //(energy/1000)^b
      par_num += 1;
      pars[par_num] = 1.0;  //c^(1000/energy)
      par_num += 1;
    }//if( input.phys_model_use_hoerl )
    
    assert( par_num == num_parameters );
    if( par_num != num_parameters )
      throw logic_error( "Num paramaters doesnt match expected" );
    
    if( !constant_parameters.empty() )
    {
      ceres::Manifold *subset_manifold = new ceres::SubsetManifold( static_cast<int>(num_parameters), constant_parameters );
      problem.SetManifold( pars, subset_manifold );
    }
  }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  
  
  //problem.AddParameterBlock( pars, static_cast<int>(num_parameters), subset_manifold ); // I guess takes ownership of subset_manifold?  TODO: check no memory leak

// Set a lower bound on relative activities to be 0
  for( size_t i = 0; i < num_nuclides; ++i )
    problem.SetParameterLowerBound( pars, static_cast<int>(i), 0.0 );
  
  if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
  {
    const auto set_bounds = [&problem, num_nuclides, pars]( const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt, size_t &index ){
      if( !opt || (!opt->material && ((opt->atomic_number < 0.999) || (opt->atomic_number > 98.001))) )
        return;

      if( opt->fit_atomic_number )
      {
        double lower_an = opt->lower_fit_atomic_number;
        double upper_an = opt->upper_fit_atomic_number;
        if( (lower_an == upper_an) && (lower_an == 0.0) )
        {
          lower_an = 1.0;
          upper_an = 98.0;
        }

        problem.SetParameterLowerBound( pars, static_cast<int>(index), lower_an / RelActCalc::ns_an_ceres_mult );
        problem.SetParameterUpperBound( pars, static_cast<int>(index), upper_an / RelActCalc::ns_an_ceres_mult );
        index += 1; //Add parameter for AN, only if fitting it
      }
      
      if( opt->fit_areal_density )
      {
        double lower_ad = opt->lower_fit_areal_density / PhysicalUnits::g_per_cm2;
        double upper_ad = opt->upper_fit_areal_density / PhysicalUnits::g_per_cm2;
        if( (lower_ad == upper_ad) && (lower_ad == 0.0) )
        {
          lower_ad = 0.0;
          upper_ad = RelActCalc::PhysicalModelShieldInput::sm_upper_allowed_areal_density_in_g_per_cm2;
        }
        
        problem.SetParameterLowerBound( pars, static_cast<int>(index), lower_ad );
        problem.SetParameterUpperBound( pars, static_cast<int>(index), upper_ad );
      }
      
      index += 1; //Add parameter for AD, always
    };//set_bounds lambda

    size_t index = num_nuclides;
    set_bounds( input.phys_model_self_atten, index );
    for( size_t i = 0; i < input.phys_model_external_attens.size(); ++i )
      set_bounds( input.phys_model_external_attens[i], index );
    
    if( input.phys_model_use_hoerl )
    {
      assert( num_parameters > 2 );
      problem.SetParameterLowerBound( pars, static_cast<int>(index), 0.0 );
      problem.SetParameterUpperBound( pars, static_cast<int>(index), 2.0 );
      index += 1;
      assert( index < num_parameters );
      problem.SetParameterLowerBound( pars, static_cast<int>(index), -3.0 );
      problem.SetParameterUpperBound( pars, static_cast<int>(index), 3.0 );
      index += 1;
    }
    
    assert( index == num_parameters );
  }else if( input.use_ceres_to_fit_eqn )
  {
    try
    {
      vector<double> rel_activities( num_nuclides ), dummy_parameters( num_parameters, 1.0 );
      for( size_t i = 0; i < num_nuclides; ++i )
        rel_activities[i] = cost_functor->relative_activity( cost_functor->m_isotopes[i], parameters );

      vector<double> fit_pars;
      fit_rel_eff_eqn_lls( eqn_form, eqn_order, cost_functor->m_isotopes, rel_activities, peak_infos, fit_pars, nullptr );
      assert( fit_pars.size() == (eqn_order + 1) );
      assert( parameters.size() == (cost_functor->m_isotopes.size() + eqn_order + 1) );
      for( size_t i = 0; i < (eqn_order + 1); ++i )
        parameters[num_nuclides + i] = fit_pars[i];
    }catch( std::exception &e )
    {
      solution.m_status = ManualSolutionStatus::ErrorInitializing;
      solution.m_error_message = "Failed to fit initial relative efficiency equation: " + string(e.what());
      return solution;
    }
  }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel ) / else if( input.use_ceres_to_fit_eqn )

  // Okay - we've set our problem up
  ceres::Solver::Options ceres_options;
  ceres_options.linear_solver_type = ceres::DENSE_QR;
  ceres_options.logging_type = ceres::SILENT;
  ceres_options.minimizer_progress_to_stdout = false; //true;
  ceres_options.max_num_iterations = 50;
  ceres_options.max_solver_time_in_seconds = 60.0;
  //ceres_options.min_trust_region_radius = 1e-10;
  
  //ceres_options.use_nonmonotonic_steps = true;
  //ceres_options.max_consecutive_nonmonotonic_steps = 5;
  // There are a lot more options that could be useful here!

  //parameter_tolerance = 1e-8;
  //double gradient_tolerance = 1e-10;
  //double function_tolerance = 1e-6;
  //int max_num_consecutive_invalid_steps = 5;
  
  // Setting ceres_options.num_threads >1 doesnt seem to do much (any?) good
  ceres_options.num_threads = std::thread::hardware_concurrency();
  if( !ceres_options.num_threads )
  {
    assert( 0 );
    ceres_options.num_threads = 4;
  }
  
  //cout << "Starting parameter values: {";
  //for( size_t i = 0; i < num_parameters; ++i )
  //  cout << (i ? ", " : "") << parameters[i];
  //cout << "}\n";
  
  ceres::Solver::Summary summary;
  
  try
  {
    ceres::Solve(ceres_options, &problem, &summary);
    
    solution.m_num_function_eval_solution = static_cast<int>( cost_functor->m_ncalls );
    
    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        // good deal, all is well with our little fit
        break;
        
      case ceres::NO_CONVERGENCE:
      case ceres::FAILURE:
      case ceres::USER_FAILURE:
        throw runtime_error( "The L-M solving failed." );
        break;
    }//switch( summary.termination_type )
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorFindingSolution;
    solution.m_error_message += e.what();
    
    cerr << "RelActCalcManual::solve_relative_efficiency: Failed in solving solution: "
         << e.what() << endl;
    
    return solution;
  }//try / catch to fit the solution
  
  
  std::cout << summary.BriefReport() << "\n";
  //std::cout << summary.FullReport() << "\n";
  const auto nmicro = std::chrono::duration<double, std::micro>(std::chrono::high_resolution_clock::now() - start_time).count();
  cout << "Took " << solution.m_num_function_eval_solution << " calls and " << setprecision(6) << nmicro << " us to solve." << endl;
  cout << "Final parameter values: {";
  for( size_t i = 0; i < num_parameters; ++i )
    cout << (i ? ", " : "") << parameters[i];
  cout << "}\n";
  cout << "Chi2=" << summary.final_cost << " (from initial value " << summary.initial_cost << ")\n\n";
  

  solution.m_fit_parameters = parameters;

  try
  {
    ceres::Covariance::Options cov_options;
    //cov_options.algorithm_type = ceres::CovarianceAlgorithmType::SPARSE_QR; //faster, but not capable of computing the covariance if the Jacobian is rank deficient.
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::DENSE_SVD;
    cov_options.null_space_rank = -1;
    cov_options.num_threads = ceres_options.num_threads;
    
    vector<double> uncertainties( num_nuclides, 0.0 ), uncerts_squared( num_nuclides, 0.0 );
    
    ceres::Covariance covariance(cov_options);
    vector<const double*> parameter_blocks;
    vector<pair<const double*, const double*> > covariance_blocks;

    parameter_blocks.push_back( pars );
    covariance_blocks.push_back( {pars, pars} );
    
    solution.m_nonlin_covariance.clear();
    if( !covariance.Compute(covariance_blocks, &problem) )
    {
      cerr << "Failed to compute final covariances!" << endl;
      solution.m_warnings.push_back( "Failed to compute final covariances." );
    }else
    {
      // row-major order: the elements of the first row are consecutively given, followed by second
      //                  row contents, etc
      vector<double> row_major_covariance( num_parameters * num_parameters );
      
      const bool success = covariance.GetCovarianceMatrix( parameter_blocks, row_major_covariance.data() );
      assert( success );
      if( !success )
        throw runtime_error( "Failed to get covariance matrix - maybe didnt add all covariance blocks?" );
      
      solution.m_nonlin_covariance.resize( num_parameters, vector<double>(num_parameters,0.0) );
      for( size_t row = 0; row < num_parameters; ++row )
      {
        for( size_t col = 0; col < num_parameters; ++col )
          solution.m_nonlin_covariance[row][col] = row_major_covariance[row*num_parameters + col];
      }//for( size_t row = 0; row < num_nuclides; ++row )
    }//if( we failed to get covariance ) / else
    
    // Compute the Jacobian
    try
    {
      const size_t num_residuals = cost_functor->number_residuals();
      vector<double> residuals( num_residuals );
      vector<double> jacobian( num_parameters * num_residuals ); 

      const double * const parameters_ptr = parameters.data();
      double * const residuals_ptr = &(residuals[0]);
      double * jacobians_ptr = &(jacobian[0]);  //We only have a single paramater block; the pointer passed into Ceres expects first index to index into parameter block

      const bool success = cost_function->Evaluate( &parameters_ptr, residuals_ptr, &jacobians_ptr );
      if( !success )
        throw std::runtime_error( "Failed to evaluate the cost function." );

      solution.m_nonlin_jacobian.resize( num_residuals, vector<double>(num_parameters, 0.0) );
      for( size_t k = 0; k < num_residuals; ++k )
      {
        for( size_t i = 0; i < num_parameters; ++i )
          solution.m_nonlin_jacobian[k][i] = jacobian[k*num_parameters + i];
      }
    }catch(const std::exception& e)
    {
      cerr << "Failed to compute final Jacobian! - " << e.what() << endl;
      solution.m_warnings.push_back( "Failed to compute Jacobian: " + string(e.what()) );
    }//try / catch to compute Jacobian


    // Compute the relative activities
    vector<double> rel_activities( num_nuclides );
    
    for( size_t i = 0; i < num_nuclides; ++i )
      rel_activities[i] = cost_functor->relative_activity( cost_functor->m_isotopes[i], parameters );
    
    if( (eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel) || input.use_ceres_to_fit_eqn )
    {
      solution.m_rel_eff_eqn_coefficients = {pars + num_nuclides, pars + num_parameters };

      if( !solution.m_nonlin_covariance.empty() )
      {
        // The covariance matrix for the relative efficiency equation is the lower-right
        //  submatrix of the full covariance matrix.
        const size_t num_rel_eff_params = num_parameters - num_nuclides;
        assert( solution.m_nonlin_covariance.size() == num_parameters );

        solution.m_rel_eff_eqn_covariance.resize( num_rel_eff_params, vector<double>(num_rel_eff_params, 0.0) );
        for( size_t i = num_nuclides; i < num_parameters; ++i )
        {
          assert( solution.m_nonlin_covariance[i].size() == num_parameters );
          for( size_t j = num_nuclides; j < num_parameters; ++j )
            solution.m_rel_eff_eqn_covariance[i-num_nuclides][j-num_nuclides] = solution.m_nonlin_covariance[i][j];
        }
      }//if( we have the covariance matrix )
    }else
    {
      fit_rel_eff_eqn_lls( eqn_form, eqn_order, cost_functor->m_isotopes, rel_activities, peak_infos,
                          solution.m_rel_eff_eqn_coefficients, &(solution.m_rel_eff_eqn_covariance) );
      assert( solution.m_rel_eff_eqn_coefficients.size() == (eqn_order + 1) );
    }
    
    for( size_t i = 0; i < cost_functor->m_isotopes.size(); ++i )
    {
      const double rel_act_norm = cost_functor->m_rel_act_norms[i];
      
      IsotopeRelativeActivity rel_act;
      rel_act.m_isotope = cost_functor->m_isotopes[i];
      
      rel_act.m_rel_activity = rel_act_norm * parameters[i];
      if( solution.m_nonlin_covariance.empty() )
      {
        rel_act.m_rel_activity_uncert = -1.0;
      }else
      {
        assert( i < solution.m_nonlin_covariance.size() );
        assert( i < solution.m_nonlin_covariance[i].size() );
        rel_act.m_rel_activity_uncert = rel_act_norm * std::sqrt( solution.m_nonlin_covariance[i][i] );
      }
      
      solution.m_rel_activities.push_back( std::move(rel_act) );
    }//for( loop over relative activities )
    
    if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      const auto get_res = [num_nuclides, &parameters, &solution]( const RelActCalc::PhysicalModelShieldInput &orig_opt, 
                                                              size_t &shield_index ) -> unique_ptr<RelEffSolution::PhysModelShieldFit> {
        if( !orig_opt.material
           && ((orig_opt.atomic_number < 1.0) || (orig_opt.atomic_number > 98.0))
           && !orig_opt.fit_atomic_number )
        {
          return nullptr;
        }
        
        assert( solution.m_nonlin_covariance.empty()
                || (solution.m_nonlin_covariance.size() > shield_index) );

        auto shield_result = std::make_unique<RelEffSolution::PhysModelShieldFit>();
        shield_result->m_material = orig_opt.material;
        if( !shield_result->m_material )
        {
          if( orig_opt.fit_atomic_number )
          {
            shield_result->m_atomic_number = parameters[shield_index] * RelActCalc::ns_an_ceres_mult;
            
            if( !solution.m_nonlin_covariance.empty()  )
            {
              assert( shield_index < solution.m_nonlin_covariance.size() );
              assert( orig_opt.fit_atomic_number
                     || (solution.m_nonlin_covariance[shield_index][shield_index] == 0.0) );
              
              if( orig_opt.fit_atomic_number )
                shield_result->m_atomic_number_uncert = sqrt( solution.m_nonlin_covariance[shield_index][shield_index] ) * RelActCalc::ns_an_ceres_mult;
            }//if( !solution.m_nonlin_covariance.empty()  )
          
            shield_index += 1;
          }else
          {
            shield_result->m_atomic_number = orig_opt.atomic_number;
          }
        }//if( !shield_result->m_material )
        
        const double ad_g_cm2 = parameters[shield_index];
        shield_result->m_areal_density = ad_g_cm2 * PhysicalUnits::g_per_cm2;

        if( !solution.m_nonlin_covariance.empty()  )
        {
          //assert( orig_opt.fit_areal_density
          //        || (solution.m_nonlin_covariance[shield_index][shield_index] == 0.0) );
          
          if( orig_opt.fit_areal_density )
            shield_result->m_areal_density_uncert = sqrt( solution.m_nonlin_covariance[shield_index][shield_index] ) * PhysicalUnits::g_per_cm2;
        }//if( !solution.m_nonlin_covariance.empty()  )
        
        shield_index += 1; //increment for areal density
        
        return std::move(shield_result);
      };//get_res lamda

      size_t shield_index = num_nuclides;
      if( input.phys_model_self_atten )
        solution.m_phys_model_self_atten_shield = get_res( *input.phys_model_self_atten, shield_index );

      for( const shared_ptr<const RelActCalc::PhysicalModelShieldInput> &opt : input.phys_model_external_attens )
      {
        assert( opt );
        if( opt )
          solution.m_phys_model_external_atten_shields.push_back( get_res( *opt, shield_index ) );
      }//for( size_t i = 0; i < input.phys_model_external_attens.size(); ++i )
      
      assert( (shield_index + (solution.m_input.phys_model_use_hoerl ? 2 : 0)) == num_parameters );
    }//if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    
    // We'll manually compute the Chi2 here, so we only take into account statistical uncertainties,
    //  (i.e., we ignore #GenericPeakInfo::m_base_rel_eff_uncert)
    solution.m_chi2 = 0.0;
    
    // TODO: The DOF is probably off by one - need to think on this and come back to it
    //assert( cost_functor->m_peak_infos.size() >= ((eqn_order+1) + (cost_functor->m_isotopes.size() - 1)) );
    if( cost_functor->m_input.peaks.size() < ((eqn_order+1) + (cost_functor->m_isotopes.size() - 1)) )
      throw runtime_error( "There are only " + std::to_string(cost_functor->m_input.peaks.size())
                          + " peaks, but you are asking to fit " + std::to_string(eqn_order+1)
                          + " rel. eff. parameters, and "
                          + std::to_string(cost_functor->m_isotopes.size())
                          + " isotope rel. act."
                          );
    
    const int num_peaks = static_cast<int>(cost_functor->m_input.peaks.size());
    const int num_isotopes = static_cast<int>(cost_functor->m_isotopes.size());
    if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
    {
      solution.m_dof = num_peaks - (num_isotopes - 1);
      if( input.phys_model_self_atten )
      {
        solution.m_dof -= static_cast<int>(input.phys_model_self_atten->fit_areal_density);
        if( !input.phys_model_self_atten->material )
          solution.m_dof -= static_cast<int>(input.phys_model_self_atten->fit_atomic_number);
      }

      for( const auto &opt : input.phys_model_external_attens )
      {
        solution.m_dof -= static_cast<int>(opt->fit_areal_density);
        if( !opt->material )
          solution.m_dof -= static_cast<int>(opt->fit_atomic_number);
      }

      solution.m_dof -= (input.phys_model_use_hoerl ? 2 : 1); // for b and c
    }else
    {
      solution.m_dof = static_cast<int>( num_peaks - (eqn_order + 1) - (num_isotopes - 1) );
    }

    
    const vector<double> &rel_eff_coefs = solution.m_rel_eff_eqn_coefficients;
    ManualGenericRelActFunctor::PhysModelRelEqnDef<double> phys_mode_rel_eqn_input;
    if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      phys_mode_rel_eqn_input = ManualGenericRelActFunctor::make_phys_eqn_input( input, rel_eff_coefs );
    
    bool used_unweighted = false, used_add_uncert = false;
    for( const GenericPeakInfo &peak : cost_functor->m_input.peaks )
    {
      double curve_val;
      if( eqn_form == RelActCalc::RelEffEqnForm::FramPhysicalModel )
      {
        curve_val = RelActCalc::eval_physical_model_eqn( peak.m_energy,
                                                  phys_mode_rel_eqn_input.self_atten,
                                                  phys_mode_rel_eqn_input.external_attens,
                                                  phys_mode_rel_eqn_input.det.get(),
                                                  phys_mode_rel_eqn_input.hoerl_b,
                                                  phys_mode_rel_eqn_input.hoerl_c );
      }else
      {
        curve_val = RelActCalc::eval_eqn( peak.m_energy, eqn_form, rel_eff_coefs );
      }
      
      used_add_uncert |= (peak.m_base_rel_eff_uncert > 1.0E-9);
      used_unweighted |= (peak.m_base_rel_eff_uncert < -1.0E-9);
      
      double expected_src_counts = 0.0;
      for( const GenericLineInfo &line : peak.m_source_gammas )
      {
        const double rel_activity = cost_functor->relative_activity( line.m_isotope, parameters );
        expected_src_counts += rel_activity * line.m_yield;
      }//for( const RelActCalc::GammaLineInfo &line : peak.m_source_gammas )
      
      const double expected_counts = expected_src_counts * curve_val;
      solution.m_chi2 += std::pow( (expected_counts - peak.m_counts) / peak.m_counts_uncert, 2.0 );
    }//for( loop over energies to evaluate at )
    
    if( used_add_uncert )
      solution.m_warnings.push_back( "Additional uncertainties were applied to peaks"
                                     " - the result uncertainties include these, so may not be"
                                     " reliable to interpret. The &chi;<sup>2</sup>/DOF does"
                                     " not include the add. uncerts." );
    
    if( used_unweighted )
      solution.m_warnings.push_back( "Fit to rel. eff. was unweighted, so uncertainties may not have much meaning." );
    
    
    solution.m_status = ManualSolutionStatus::Success;
    solution.m_num_function_eval_total = static_cast<int>( cost_functor->m_ncalls );
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorGettingSolution;
    solution.m_error_message = e.what();
    cerr << "RelActCalcManual::solve_relative_efficiency: Failed to get solution after solving: "
         << e.what() << endl;
    
    return solution;
  }//try / catch to get the solution
  
  /*
  if( (solution.m_status == ManualSolutionStatus::Success)
     && !solution.m_nonlin_covariance.empty() && !solution.m_fit_parameters.empty() )
  {
    cout << "Covariance matrix:" << endl;
    cout << "  " << setw(10) << " " << " ";
    for( size_t row = 0; row < num_parameters; ++row )
      cout << setw(10) << solution.parameter_name(row) << " ";
    cout << endl;
    
    for( size_t row = 0; row < num_parameters; ++row )
    {
      cout << "  " << setw(10) << solution.parameter_name(row) << " ";
      for( size_t col = 0; col < num_parameters; ++col )
        cout << setw(10) << setprecision(2) << solution.m_nonlin_covariance[row][col] << " ";
      cout << endl;
    }//for( size_t row = 0; row < num_nuclides; ++row )
  }//if( we have the covariance matrix )
  */
  
  return solution;
}//solve_relative_efficiency(...)




namespace PeakCsvInput
{
NuclideInfo::NuclideInfo( const char *p, const char *nuc, bool opt, float kev, float br )
: parent(p), source_nuclide(nuc), energy(kev), yield(br), optional(opt)
{
}

const char *to_str( const NucDataSrc src )
{
  switch( src )
  {
    case NucDataSrc::Icrp107_U:        return "Icrp107_U";
    case NucDataSrc::Lanl_U:           return "Lanl_U";
    case NucDataSrc::IcrpLanlGadras_U: return "IcrpLanlGadras_U";
    case NucDataSrc::SandiaDecay:      return "SandiaDecay";
    case NucDataSrc::Undefined:        return "Undefined";
  }
  assert( 0 );
  return "";
}//to_str( NucDataSrc )



NucMatchResults fill_in_nuclide_info( const vector<RelActCalcManual::GenericPeakInfo> peaks,
                                     const NucDataSrc nuc_data_src,
                                     const vector<pair<float,float>> energy_ranges,
                                     std::vector<NucAndAge> isotopes,
                                     const float energy_tolerance_sigma,
                                     const vector<float> excluded_peak_energies,
                                     const float measurement_duration )
{
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  assert( db );
  
  //https://journals.sagepub.com/doi/pdf/10.1177/ANIB_38_3
  // and more specifically the supplementary material at
  // https://journals.sagepub.com/doi/suppl/10.1177/ANIB_38_3
  const vector<NuclideInfo> icrp107{{//18 entries
    {"U235", "U235",         false, 143.76f,      0.11f},
    {"U235", "U235",         false, 163.33f,      0.0508f},
    {"U235", "U235",         false, 185.715f,     0.572f},
    {"U235", "U235",         true,  202.11f,      0.0108f},
    {"U235", "U235",         false, 205.311f,     0.0501f},
    
    {"U232", "Pb212",        false, 238.632f,     0.433f},
    {"U232", "Tl208",        false, 583.191f,     0.304f},
    {"U232", "Bi212",        false, 727.33f,      0.0658f},
    {"U232", "Tl208",        false, 860.564f,     0.0447f},
    
    {"U238", "Pa234m",       false, 258.26f,      0.000726833f},
    {"U238", "Pa234",        true,  569.3173913f, 0.00018952f},
    {"U238", "Pa234/Pa234m", true,  742.81f,      0.000831678f},
    {"U238", "Pa234m",       false, 766.36f,      0.002935286f},
    {"U238", "Pa234/Pa234m", true,  880.5742495f, 0.000204387f},
    {"U238", "Pa234/Pa234m", true,  883.24f,      0.000188208f},
    {"U238", "Pa234/Pa234m", true,  945.9684295f, 0.000313081f},
    {"U238", "Pa234m",       false, 1001.03f,     0.008356588f},
    
    {"U234", "U234",         false, 120.9f,       0.000397362f},
  }};//icrp107
  
  
  const vector<NuclideInfo> lanl{{ //17 entries
    {"U235", "U235",          false, 143.76f,  0.11f},
    {"U235", "U235",          false, 163.36f,  0.0505f},
    {"U235", "U235",          false, 185.715f, 0.57f},
    {"U235", "U235",          true,  202.11f,  0.01098f},
    {"U235", "U235",          false, 205.311f, 0.0503f},
    
    {"U232", "Pb212",         false, 238.625f, 0.48f},
    {"U232", "Tl208",         false, 583.187f, 0.306f},
    {"U232", "Bi212",         false, 727.3f,   0.0676f},
    {"U232", "Tl208",         false, 860.56f,  0.046f},
    
    {"U238", "Pa234m",        false, 258.26f,  0.000754f},
    {"U238", "Pa234/Pa234m",  true,  742.83f,  0.000907f},
    {"U238", "Pa234m",        false, 766.4f,   0.003074f},
    {"U238", "Pa234/Pa234m",  true,  880.47f,  0.000213f},
    {"U238", "Pa234/Pa234m",  true,  883.24f,  0.000213f},
    {"U238", "Pa234/Pa234m",  true,  945.95f,  0.000347f},
    {"U238", "Pa234m",        false, 1001.03f, 0.008371f},
    
    {"U234", "U234",          false, 120.905f, 0.00035f}
  }};
  
  
  const vector<NuclideInfo> lanl_icrp_gad{{ //21 entries
    {"U235", "U235",          false, 143.76f,      0.11f},
    {"U235", "U235",          false, 163.36f,      0.0505f},
    {"U235", "U235",          false, 185.715f,     0.57f},
    {"U235", "U235",          true,  202.11f,      0.01098f},
    {"U235", "U235",          false, 205.311f,     0.0503f},
    {"U235", "U235",          true,  221.38f,      0.0012f},
    {"U235", "U235",          true,  246.84f,      0.00053f},
    {"U235", "U235",          true,  345.9f,       0.0003f},
    
    {"U232", "Pb212",         false, 238.625f,     0.48f},
    {"U232", "Tl208",         false, 583.187f,     0.306f},
    {"U232", "Bi212",         false, 727.3f,       0.0676f},
    {"U232", "Tl208",         false, 860.56f,      0.046f},
    
    {"U238", "Pa234m",        false, 258.26f,      0.000754f},
    {"U238", "Pa234",         true,  569.3173913f, 0.00018952f},
    {"U238", "Pa234/Pa-234m", true,  742.83f,      0.000907f},
    {"U238", "Pa234m",        false, 766.4f,       0.003074f},
    {"U238", "Pa234/Pa-234m", true,  880.47f,      0.000213f},
    {"U238", "Pa234/Pa-234m", true,  883.24f,      0.000213f},
    {"U238", "Pa234/Pa-234m", true,  945.95f,      0.000347f},
    {"U238", "Pa234m",        false, 1001.03f,     0.008371f},
    
    {"U234", "U234",          false, 120.905f,     0.00035f}
  }};
  
  // Check isotopes are valid, and normalize their name, and if there is an age, get it.
  if( isotopes.empty() )
  {
    if( (nuc_data_src == NucDataSrc::SandiaDecay) || (nuc_data_src == NucDataSrc::Undefined) )
      throw runtime_error( "fill_in_nuclide_info: No nuclides, or specialized nuclear data sources specified." );
       
    vector<NuclideInfo> specialize_src;
    switch( nuc_data_src )
    {
      case NucDataSrc::Icrp107_U:
        specialize_src = icrp107;
        break;
        
      case NucDataSrc::Lanl_U:
        specialize_src = lanl;
        break;
        
      case NucDataSrc::IcrpLanlGadras_U:
        specialize_src = lanl_icrp_gad;
        break;
        
      case NucDataSrc::SandiaDecay:
      case NucDataSrc::Undefined:
        assert( 0 );
        break;
    }//switch( nuc_data_src )
    
    
    for( const auto &info : specialize_src )
    {
      bool has_iso = false;
      for( const auto &iso : isotopes )
        has_iso |= (iso.nuclide == info.parent);
      
      if( !has_iso )
      {
        const SandiaDecay::Nuclide *nuc = db->nuclide(info.parent);
        assert( nuc );
        if( !nuc )
          throw runtime_error( "Some how '" + info.parent + "' isnt a valid nuclide." );
        
        isotopes.emplace_back( nuc->symbol, -1.0, false );
      }//
    }//for( const auto &info : specialize_src )
  }//if( isotopes.empty() )
  
  
  // We will put raw source info into 'initial_nucs_info', and then filter this down to the
  //  nuclides and energy ranges we will actually use
  //  If we are correcting activities for nuclides decay during the measurement, then we will
  //  place the uncorrected info into `non_decay_corr_initial_nucs_info` (that is,
  //  `non_decay_corr_initial_nucs_info` will not have any entries if we arent decay correcting)
  vector<NuclideInfo> initial_nucs_info, non_decay_corr_initial_nucs_info;
  for( size_t i = 0; i < isotopes.size(); ++i )
  {
    string &iso = isotopes[i].nuclide;
    double &age = isotopes[i].age;
    
    const ReactionGamma::Reaction *rctn = nullptr;
    const SandiaDecay::Nuclide *nuc = db->nuclide( iso );
    
    if( !nuc )
    {
      const ReactionGamma *reactiondb = ReactionGammaServer::database();
      
      try
      {
        vector<ReactionGamma::ReactionPhotopeak> reactions;
        reactiondb->gammas( iso, reactions );
        if( reactions.empty() )
          throw runtime_error( "unknown reaction" );
        
        // TODO: currently just using first possible reaction; need to implement retrieving reactions by name from ReactionGamma.
        rctn = reactions[0].reaction;
      }catch( std::exception &e )
      {
        cerr << "Invalid reaction ("<< iso << "): " << e.what() << endl;
      }
    }//if( !nuc )
    
    if( !nuc && !rctn )
      throw runtime_error( "Invalid nuclide '" + iso + "' specified." );
    
    iso = nuc ? nuc->symbol : rctn->name();
    
    // Make sure we try to use a U-data source, only for U, and the nuclides of it we have,
    //  otherwise we'll default back to using SandiaDecay
    auto src = nuc_data_src;
    if( (nuc && (nuc->atomicNumber != 92)) || rctn )
    {
      src = NucDataSrc::SandiaDecay;
    }else if( nuc->atomicNumber == 92 )
    {
      auto has_U_nuc = [db]( const SandiaDecay::Nuclide *nuc, const vector<NuclideInfo> &input ) -> bool {
        for( const NuclideInfo &inputnuc : input )
        {
          const SandiaDecay::Nuclide *inputnuc_ptr = db->nuclide( inputnuc.parent );
          assert( inputnuc_ptr );
          if( inputnuc_ptr == nuc )
            return true;
        }//for( const NuclideInfo &inputnuc : input )
        return false;
      };//filter_for_nuc(...)
      
      // TODO: we should probably give a warning that we are using SandiaDecay data if the
      //       uranium dataset doesn't have this nuclide.
      switch( src )
      {
        case NucDataSrc::Icrp107_U:
          if( !has_U_nuc(nuc, icrp107) )
            src = NucDataSrc::SandiaDecay;
        break;
          
        case NucDataSrc::Lanl_U:
          if( !has_U_nuc(nuc, lanl) )
            src = NucDataSrc::SandiaDecay;
        break;
          
        case NucDataSrc::IcrpLanlGadras_U:
          if( !has_U_nuc(nuc, lanl_icrp_gad) )
            src = NucDataSrc::SandiaDecay;
        break;
          
        case NucDataSrc::SandiaDecay:
        case NucDataSrc::Undefined:
          break;
      }//switch( src )
    }//if( non-U nuclide ) / else
    
    
    switch( src )
    {
      case NucDataSrc::Icrp107_U:
      case NucDataSrc::Lanl_U:
      case NucDataSrc::IcrpLanlGadras_U:
      case NucDataSrc::Undefined:
        age = -1.0;
        break;
        
      case NucDataSrc::SandiaDecay:
        if( nuc && (age < 0.0) )
          age = PeakDef::defaultDecayTime(nuc, nullptr);
        break;
    }//switch( nucdatasrc )
    
    
    auto filter_for_nuc = [db]( const SandiaDecay::Nuclide *nuc, const vector<NuclideInfo> &input ) -> vector<NuclideInfo> {
      vector<NuclideInfo> results;
      
      for( const NuclideInfo &inputnuc : input )
      {
        const SandiaDecay::Nuclide *inputnuc_ptr = db->nuclide( inputnuc.parent );
        assert( inputnuc_ptr );
        if( inputnuc_ptr == nuc )
          results.push_back( inputnuc );
      }//for( const NuclideInfo &inputnuc : input )
      
      return results;
    };//filter_for_nuc(...)
    
    
    switch( src )
    {
      case NucDataSrc::Icrp107_U:
      {
        const auto these_infos = filter_for_nuc( nuc, icrp107 );
        initial_nucs_info.insert( end(initial_nucs_info), begin(these_infos), end(these_infos) );
        break;
      }
        
      case NucDataSrc::Lanl_U:
      {
        const auto these_infos = filter_for_nuc( nuc, lanl );
        initial_nucs_info.insert( end(initial_nucs_info), begin(these_infos), end(these_infos) );
        break;
      }
        
      case NucDataSrc::IcrpLanlGadras_U:
      {
        const auto these_infos = filter_for_nuc( nuc, lanl_icrp_gad );
        initial_nucs_info.insert( end(initial_nucs_info), begin(these_infos), end(these_infos) );
        break;
      }
        
      case NucDataSrc::SandiaDecay:
      {
        if( isotopes.empty() )
          throw runtime_error( "You must specify the isotopes to use when using SandiaDecay as your source data." );
        
        if( rctn )
        {
          for( const auto &g : rctn->gammas )
          {
            if( g.abundance > std::numeric_limits<float>::min() )
              initial_nucs_info.emplace_back( rctn->name().c_str(), "", true, g.energy, g.abundance );
          }
        }else
        {
          assert( nuc );
          const double ref_act = 1.0*SandiaDecay::MBq;
          const double decrease_factor = std::exp( -age * nuc->decayConstant() );
          const double initial_activity = ref_act / decrease_factor;
          SandiaDecay::NuclideMixture mix;
          mix.addNuclideByActivity( nuc, initial_activity );
          
          assert( fabs(ref_act - mix.activity(age, nuc)) < 0.001*ref_act );
          
          vector<SandiaDecay::EnergyRatePair> un_decay_corrected_photons;
          vector<SandiaDecay::EnergyRatePair> photons = mix.photons( age, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy );
          
          if( isotopes[i].decay_during_measurement )
          {
            vector<SandiaDecay::EnergyRatePair> corr_photons
                  = GammaInteractionCalc::decay_during_meas_corrected_gammas( mix, age, measurement_duration );
            
            // See what the decay correction is here, and note to then include in the results.
            assert( corr_photons.size() == photons.size() );
            assert( std::is_sorted(begin(corr_photons), end(corr_photons),
                   []( const SandiaDecay::EnergyRatePair &lhs, const SandiaDecay::EnergyRatePair &rhs ){
              return lhs.energy < rhs.energy;
            }) );
            assert( std::is_sorted(begin(photons), end(photons),
                   []( const SandiaDecay::EnergyRatePair &lhs, const SandiaDecay::EnergyRatePair &rhs ){
              return lhs.energy < rhs.energy;
            }) );
            //for( size_t i = 0; i < std::min(photons.size(),corr_photons.size()); ++i ){
            //  assert( fabs(photons[i].energy - corr_photons[i].energy) < 0.001 );
            //}
            
            un_decay_corrected_photons = std::move(photons);
            photons = std::move(corr_photons);
          }//if( isotopes[i].decay_during_measurement ) / else
          
          for( size_t info_index = 0; info_index < photons.size(); ++info_index )
          {
            const SandiaDecay::EnergyRatePair &rate_info = photons[info_index];
              
            const char * const parent = nuc->symbol.c_str();
            const char * const resp_nuc = ""; //TODO: bother to get the nuclide thats actually giving off this gamma; see decay_gammas(...) in RelActCalcAuto.cpp
            const bool optional_use = true;
            const float energy = rate_info.energy;
            const float br = static_cast<float>( rate_info.numPerSecond / ref_act );
            
            if( br > std::numeric_limits<float>::min() )
            {
              initial_nucs_info.emplace_back( parent, resp_nuc, optional_use, energy, br );
              
              if( info_index < un_decay_corrected_photons.size() )
              {
                assert( un_decay_corrected_photons.size() == photons.size() );
                const SandiaDecay::EnergyRatePair &noncorr_rate_info = un_decay_corrected_photons[info_index];
                assert( fabs(noncorr_rate_info.energy - energy) < 0.001 );
                
                const float non_corr_br = static_cast<float>( noncorr_rate_info.numPerSecond / ref_act );
                non_decay_corr_initial_nucs_info.emplace_back( parent, resp_nuc, optional_use, energy, non_corr_br );
              }//if( decay correct nuclides )
            }//if( br > std::numeric_limits<float>::min() )
          }//for( const SandiaDecay::EnergyRatePair info : photons )
        }//if( rctn ) / else( nuc )
        
        break;
      }//case NucDataSrc::SandiaDecay:
        
      case NucDataSrc::Undefined:
        assert(0);
        throw std::logic_error( "shouldnt be here" );
        break;
    }//switch( nucdatasrc )
  }//for( size_t i = 0; i < isotopes.size(); ++i )
  
  
  std::sort( begin(initial_nucs_info), end(initial_nucs_info), []( const NuclideInfo &lhs, const NuclideInfo &rhs ){
    return lhs.energy < rhs.energy;
  });

  std::sort( begin(non_decay_corr_initial_nucs_info), end(non_decay_corr_initial_nucs_info),
            []( const NuclideInfo &lhs, const NuclideInfo &rhs ){
    return lhs.energy < rhs.energy;
  });
  
  
  // Now we'll put the entries from 'initial_nucs_info' we *might* actually use, into
  //  'candidate_nucs_info'.  Note that we're doing this in stages to provide warnings to the
  //  user about peaks not used, or matched, or whatever.
  vector<NuclideInfo> candidate_nucs_info;
  for( const auto &info : initial_nucs_info )
  {
    // Check if isotope is wanted
    bool nuc_wanted = false;
    for( const auto &n : isotopes )
      nuc_wanted |= (n.nuclide == info.parent);
    
    if( !nuc_wanted )
      continue;
    
    // Check if in energy range  (energy_ranges)
    bool in_energy_range = energy_ranges.empty();
    for( const pair<float,float> &r : energy_ranges )
      in_energy_range |= ((info.energy >= r.first) && (info.energy <= r.second));
    
    if( !in_energy_range )
      continue;
    
    candidate_nucs_info.push_back( info );
  }//for( const auto &info : initial_nucs_info )
  
  // Now go through and actually match up the info
  vector<RelActCalcManual::GenericPeakInfo> matched_peaks = peaks;
  vector<bool> used_isotope( isotopes.size(), false );
  vector<bool> used_peak( matched_peaks.size(), false );
  vector<bool> peak_was_excluded( matched_peaks.size(), false );
  vector<bool> used_candidate_nucs_info( candidate_nucs_info.size(), false );
  
  // Lets keep track of peaks with gammas that have received decay-during-measurement
  //  corrections; following variable will only have un-corrected values for peaks that
  //  received corrections and only entries in #GenericPeakInfo::m_source_gammas that
  //  have been corrected.
  vector<RelActCalcManual::GenericPeakInfo> un_corrected_peaks;
  
  for( size_t peak_index = 0; peak_index < matched_peaks.size(); ++peak_index )
  {
    RelActCalcManual::GenericPeakInfo &peak = matched_peaks[peak_index];
    
    const double peak_sigma = (peak.m_fwhm > 0.0) ? (peak.m_fwhm / 2.35482) : 1.0;
    
    // Check if this peak is specifically excluded via the 'exclude-peak' command line argument
    bool exclude = false;
    for( size_t i = 0; !exclude && (i < excluded_peak_energies.size()); ++i )
      exclude = (fabs(excluded_peak_energies[i] - peak.m_energy) <= (peak_sigma * energy_tolerance_sigma) );
    peak_was_excluded[peak_index] = exclude;
    if( exclude )
      continue;
    
    RelActCalcManual::GenericPeakInfo un_decay_corr_peak = peak;
    assert( un_decay_corr_peak.m_source_gammas.empty() );
    un_decay_corr_peak.m_source_gammas.clear();
    
    // Try to match this peak to a source gamma line
    for( size_t nuc_index = 0; nuc_index < candidate_nucs_info.size(); ++nuc_index )
    {
      NuclideInfo &nuc = candidate_nucs_info[nuc_index];
      
      const double peak_sigma = (peak.m_fwhm > 0.0) ? (peak.m_fwhm / 2.35482) : 1.0;
      
      if( fabs(nuc.energy - peak.m_energy) <= (peak_sigma * energy_tolerance_sigma) )
      {
        used_peak[peak_index] = true;
        used_candidate_nucs_info[nuc_index] = true;
        
        size_t iso_pos = isotopes.size();
        for( size_t i = 0; i < isotopes.size(); ++i )
        {
          if( isotopes[i].nuclide == nuc.parent )
          {
            iso_pos = i;
            break;
          }
        }

        assert( iso_pos != isotopes.size() );
        if( iso_pos == isotopes.size() )//Shouldnt ever happen
          throw std::logic_error( "Failed to find source nuclide in isotopes to use, after filtering" );
        
        used_isotope[iso_pos] = true;
        
        bool nuc_already_used = false;
        for( RelActCalcManual::GenericLineInfo &gamma : peak.m_source_gammas )
        {
          if( gamma.m_isotope == nuc.parent )
          {
            nuc_already_used = true;
            gamma.m_yield += nuc.yield;
            break;
          }
        }//for( loop over existing source gammas for `peak` )
        
        if( !nuc_already_used )
          peak.m_source_gammas.emplace_back( nuc.yield, nuc.parent );
        
        
        //look through `non_decay_corr_initial_nucs_info` and try to match to `isotopes[i]`
        //  TODO: `non_decay_corr_initial_nucs_info` is sorted by energy, so we could be much
        //         smarter here
        for( const NuclideInfo &un_corr_nuc : non_decay_corr_initial_nucs_info )
        {
          if( (nuc.parent == un_corr_nuc.parent)
             && (nuc.source_nuclide == un_corr_nuc.source_nuclide)
             && (fabs(nuc.energy - un_corr_nuc.energy) < 0.001) )
          {
            bool uncorr_nuc_already_used = false;
            for( RelActCalcManual::GenericLineInfo &gamma : un_decay_corr_peak.m_source_gammas )
            {
              if( gamma.m_isotope == nuc.parent )
              {
                uncorr_nuc_already_used = true;
                gamma.m_yield += un_corr_nuc.yield;
                break;
              }
            }//for( loop over existing source gammas for `peak` )
            
            if( !uncorr_nuc_already_used )
              un_decay_corr_peak.m_source_gammas.emplace_back( un_corr_nuc.yield, un_corr_nuc.parent );
          }//if( we found un-corrected line corresponding to `nuc` )
          
          if( un_corr_nuc.energy > nuc.energy )
            break;
        }//for( loop over non_decay_corr_initial_nucs_info )
      }//if( peak matched source data within tolerance )
    }//for( size_t nuc_index = 0; nuc_index < candidate_nucs_info.size(); ++candidate_nucs_info )
    
    if( !un_decay_corr_peak.m_source_gammas.empty() )
      un_corrected_peaks.push_back( un_decay_corr_peak );
  }//for( size_t peak_index = 0; peak_index < matched_peaks.size(); ++peak_index )
  
  
  NucMatchResults results;
  results.data_source = nuc_data_src;
  
  results.match_sigma_tolerance = energy_tolerance_sigma;
  results.energy_ranges = energy_ranges;
  results.not_decay_corrected_peaks = un_corrected_peaks;
  
  for( size_t peak_index = 0; peak_index < matched_peaks.size(); ++peak_index )
  {
    if( peak_was_excluded[peak_index] )
    {
      results.peaks_excluded.push_back( matched_peaks[peak_index] );
    }else
    {
      auto &place = used_peak[peak_index] ? results.peaks_matched : results.peaks_not_matched;
      place.push_back( matched_peaks[peak_index] );
    }
  }//for( size_t peak_index = 0; peak_index < matched_peaks.size(); ++peak_index )
  
  
  assert( used_isotope.size() == isotopes.size() );
  for( size_t index = 0; index < used_isotope.size(); ++index )
  {
    if( used_isotope[index] )
    {
      results.used_isotopes.push_back( isotopes[index].nuclide );
      results.used_isotope_ages.push_back( isotopes[index].age );
    }else
    {
      results.unused_isotopes.push_back( isotopes[index].nuclide );
    }
  }//for( size_t index = 0; index < used_isotope.size(); ++index )
  
  
  assert( used_candidate_nucs_info.size() == candidate_nucs_info.size() );
  for( size_t index = 0; index < candidate_nucs_info.size(); ++index )
  {
    auto &place = used_candidate_nucs_info[index] ? results.source_gammas_used : results.source_gammas_not_used;
    place.push_back( candidate_nucs_info[index] );
  }//for( size_t index = 0; index < candidate_nucs_info.size(); ++index )
  
  
  return results;
}//fill_in_nuclide_info(...)



vector<RelActCalcManual::GenericPeakInfo> peak_csv_to_peaks( istream &csv )
{
  //Gadras: "Energy(keV),sigma,Rate(cps),sigma,FWHM(keV),sigma,Leakage(1/s),Centroid,FileName,RecordIdx,Title,DateTime"
  //PeakEasy: "Centroid,  Net_Area,   Net_Area,      Peak, FWHM,   FWHM,Reduced, ROI_Total,ROI"
  //          "keV,    Counts,Uncertainty,       CPS,  keV,Percent,Chi_Sqr,    Counts,ID#,  File, LiveTime, Date, Time"
  
  enum class PeakCsvFormat{ PeakEasy, Gadras, Unknown };
  
  // If the CSV file is from InterSpec, and there are any nuclides provided, we will use them.
  bool contained_nuc_ids = false;
  PeakCsvFormat csv_format = PeakCsvFormat::Unknown;
  
  
  auto split_and_trim = []( vector<string> &fields, const string &line ){
    SpecUtils::split_no_delim_compress( fields, line, "," );
    for( string &s : fields )
      SpecUtils::trim( s );
  };
  
  string line;
  while( std::getline( csv, line ) )
  {
    SpecUtils::trim( line );
    if( line.empty() || line[0] == '#' )
      continue;
    
    // Line should either be "Energy(keV),sigma,Rate(cps)...", or "Centroid,  Net_Area,   Net_Area"
    vector<string> fields;
    split_and_trim( fields, line );
    if( fields.empty() )
      continue;
    
    if( fields.size() < 9 )
      throw runtime_error( "Invalid Peak CSV header line: '" + line + "'" );
    
    if( (fields[0] == "Energy(keV)") && (fields[1] == "sigma")
       && (fields[2] == "Rate(cps)") && (fields[3] == "sigma") )
    {
      csv_format = PeakCsvFormat::Gadras;
      break;
    }else if( (fields[0] == "Centroid") && (fields[1] == "Net_Area")
             && (fields[2] == "Net_Area") && (fields[3] == "Peak") )
    {
      if( !std::getline( csv, line ) )
        throw runtime_error( "Failed to get second line of PeakEasy CSV" );
      
      split_and_trim( fields, line );
      
      if( (fields.size() < 9) || (fields[0] != "keV") || (fields[1] != "Counts")
         || (fields[2] != "Uncertainty") || (fields[3] != "CPS") )
      {
        throw runtime_error( "Second line of PeakEasy CSV file is not correct: '" + line + "'" );
      }
      
      csv_format = PeakCsvFormat::PeakEasy;
      break;
    }else
    {
      throw runtime_error( "Invalid peak CSV line: '" + line + "'" );
    }
  }//while( std::getline( csv, line ) )
  
  
  size_t energy_index, amplitude_index, amplitude_sigma_index, nuclide_index, nuclide_energy_index, fwhm_kev_index;
  
  switch( csv_format )
  {
    case PeakCsvFormat::PeakEasy:
      energy_index = 0;
      amplitude_index= 1;
      amplitude_sigma_index = 2;
      nuclide_index = 13;
      nuclide_energy_index = 14;
      fwhm_kev_index = 4;
      //const size_t peak_cps_index = 3, fwhm_percent_index = 5;
      //const size_t roi_total_counts_index = 6, roi_id_index = 7, filename_index = 8;
      //const size_t live_time_index = 9, date_index = 10, time_index = 11;
      break;
      
    case PeakCsvFormat::Gadras:
      energy_index = 0;
      amplitude_index = 2;
      amplitude_sigma_index = 3;
      nuclide_index = nuclide_energy_index = 0;
      
      //const size_t sigma_index = 1; //uncert in energy
      fwhm_kev_index = 4;//, fwhm_sigma_index = 5, leakage_per_second_index = 6;
      //const size_t centroid_index = 7, filename_index = 8, record_idx_index = 9;
      //const size_t title_index = 10, date_time_index = 11;
      break;
      
    case PeakCsvFormat::Unknown:
      throw runtime_error( "Not a peak CSV file." );
  }//switch( csv_format )
  
  
  vector<RelActCalcManual::GenericPeakInfo> answer;
  
  while( std::getline(csv, line) )
  {
    SpecUtils::trim(line);
    if( line.empty() || line[0]=='#' || (!isdigit(line[0]) && line[0]!='+' && line[0]!='-') )
      continue;
    
    vector<string> fields;
    split_and_trim( fields, line );
    
    const size_t nfields = fields.size();
    if( nfields == 0 )
      continue;
    
    if( nfields < 9 )
      throw runtime_error( "Encountered line in GADRAS CSV file with only "
                          + std::to_string(nfields) + " fields.\n\tLine: \"" + line + "\"" );
    
    try
    {
      RelActCalcManual::GenericPeakInfo info;
      info.m_mean = info.m_energy = std::stod( fields[energy_index] );
      info.m_fwhm = stod( fields[fwhm_kev_index] );
      info.m_counts = std::stod( fields[amplitude_index] );
      info.m_counts_uncert = std::stod( fields[amplitude_sigma_index] );
      
      //cout << "Found peak at " << info.m_energy << " keV"
      //     << ", with Amp=" << info.m_counts << " +- " << info.m_counts_uncert
      //     << endl;
      
      if( nuclide_index
         && nuclide_energy_index
         && (fields.size() > nuclide_energy_index)
         && !fields[nuclide_index].empty()
         && !fields[nuclide_energy_index].empty() )
      {
        contained_nuc_ids = true;
        
        const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
        
        const SandiaDecay::Nuclide *nuc = db->nuclide( fields[nuclide_index] );
        if( !nuc )
          throw runtime_error( "Invalid Nuclide ID: " + fields[nuclide_index] );
        
        const double nuc_energy = stod( fields[nuclide_energy_index] );
        info.m_energy = nuc_energy;
        
        const double age = PeakDef::defaultDecayTime( nuc, nullptr );
        
        const double fwhm = info.m_fwhm;
        
        //size_t transition_index = 0;
        //const SandiaDecay::Transition *transition = nullptr;
        //PeakDef::SourceGammaType sourceGammaType;
        //PeakDef::findNearestPhotopeak( nuc, nuc_energy, -1.0, false,
        //                               transition, transition_index, sourceGammaType );
        
        SandiaDecay::NuclideMixture mix;
        mix.addAgedNuclideByActivity( nuc, 1.0, age);
        const vector<SandiaDecay::EnergyRatePair> gammas = mix.gammas(0, SandiaDecay::NuclideMixture::HowToOrder::OrderByEnergy, false);
        double yield = 0.0;
        for( const auto &erp : gammas )
        {
          // We'll sum all gammas within 1 FWHM of the specified gamma energy; not perfect, but good enough for now.
          if( fabs(erp.energy - nuc_energy) < fwhm )
            yield += erp.numPerSecond;
        }
        
        info.m_source_gammas.push_back( {yield, nuc->symbol} );
      }//if( we have nuclide ID from file )
      
      answer.push_back( info );
    }catch( std::exception &e )
    {
      throw runtime_error( "Invalid value on line '" + line + "', " + string(e.what()) );
    }//try / catch to parse a line into a peak
  }//while( SpecUtils::safe_get_line(csv, line, 2048) )
  
  if( answer.empty() )
    throw runtime_error( "No peak rows found in file." );
  
  return answer;
}//vector<GenericPeakInfo> peak_csv_to_peaks( istream &csv )

}//namespace PeakCsvInput

}//namespace RelActCalcManual
