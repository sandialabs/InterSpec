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
#include <iostream>

#include "Eigen/Dense"
#include "ceres/ceres.h"

#include "SandiaDecay/SandiaDecay.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalc.h"
#include "InterSpec/RelActCalcManual.h"
#include "InterSpec/DecayDataBaseServer.h"
#include "InterSpec/GammaInteractionCalc.h"

using namespace std;

namespace
{
struct DoWorkOnDestruct
{
  std::function<void()> m_worker;
  DoWorkOnDestruct( std::function<void()> &&worker ) : m_worker( std::move(worker) ){}
  ~DoWorkOnDestruct(){ if(m_worker) m_worker(); }
};//struct DoWorkOnDestruct


#define USE_RESIDUAL_TO_BREAK_DEGENERACY 0

/** Functor for minimizing the relative activities; relative efficiency is fit for each set of
 activities via the #fit_rel_eff_eqn_lls function.
 */
struct ManualGenericRelActFunctor  /* : ROOT::Minuit2::FCNBase() */
{
  /** The form of relative efficiency equation to use. */
  const RelActCalc::RelEffEqnForm m_eqn_form;
  
  /** The order of relative efficiency equation to use (not equation will have one more than this
   value coefficients)
   */
  const int m_eqn_order;
  
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
  std::vector<RelActCalcManual::GenericPeakInfo> m_peak_infos;
  
  /** Warnings from setting up the problem - does not include problems evaluating things. */
  std::vector<std::string> m_setup_warnings;
  
  /** just for debug purposes, we'll keep track of how many times the eval function gets called. */
  mutable std::atomic<size_t> m_ncalls;
  
  /** Constructor for this functior.
   
   Will throw exception on error.
   */
  ManualGenericRelActFunctor( const RelActCalc::RelEffEqnForm eqn_form,
                    const int eqn_order,
                    const std::vector<RelActCalcManual::GenericPeakInfo> &peak_infos )
  : m_eqn_form( eqn_form ),
  m_eqn_order( eqn_order ),
  m_peak_infos( peak_infos ),
  m_ncalls( 0 )
  {
    if( peak_infos.size() < 2 )
      throw runtime_error( "ManualGenericRelActFunctor: you must use at least two peaks." );
    
    // Apply some sanity checks to the eqn_order.  Realistically eqn_order should probably be
    //  between 3 and 6, but we'll allow an arbitrary amount of slop here.
    if( (eqn_order < 0) || (eqn_order >= 10) )
      throw runtime_error( "ManualGenericRelActFunctor: equation order must be at least 1 and less than 10." );
    
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
    std::sort( begin(m_peak_infos), end(m_peak_infos),
              []( const RelActCalcManual::GenericPeakInfo &lhs, const RelActCalcManual::GenericPeakInfo &rhs ) -> bool {
      return lhs.m_energy < rhs.m_energy;
    } );
#endif
    
    const size_t num_peaks = peak_infos.size();
    
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
    for( size_t peak_index = 0; peak_index < peak_infos.size(); ++peak_index )
    {
      const RelActCalcManual::GenericPeakInfo &peak = peak_infos[peak_index];
      
      if( peak.m_source_gammas.empty() )
        throw std::runtime_error( "ManualGenericRelActFunctor: Peak at " + std::to_string(peak.m_energy)
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
        
        // 32 bit floats have ~7 significant figures, which is probably more than we know any
        //  branching ratios to, so we'll use a floats precisions at 1.0 (i.e. FLT_EPSILON which
        //  is the difference between 1 and the smallest floating point number greater than 1) as
        //  the minimum yield we'll consider, otherwise you shouldnt be saying this isotope
        //  contributes to the peak.
        //  Although note that not allowing the yield to be zero is arbitrary - the stuff would work
        //  fine if we let zero through.
        //  TODO: look at the various decay chains to determine better minimal value (~10E-14?).
        const double min_allowable_yield = std::numeric_limits<float>::epsilon();
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
  
    vector<double> rel_act_norm_uncerts;
    const vector<double> flat_rel_eff_coefs{ (eqn_form==RelActCalc::RelEffEqnForm::LnX ? 1.0 : 0.0), 0.0 };
    RelActCalcManual::fit_act_to_rel_eff( eqn_form, flat_rel_eff_coefs,
                       m_isotopes, m_peak_infos,
                       m_rel_act_norms, rel_act_norm_uncerts );
    
    //cout << "For initial rel act:\n";
    //for( size_t peak_index = 0; peak_index < m_peak_infos.size(); ++peak_index )
    //{
    //  const auto &peak = m_peak_infos[peak_index];
    //  assert( peak.m_source_gammas.size() == 1 ); //just for testing
    //  const RelActCalcManual::GenericLineInfo &line = peak.m_source_gammas[0];
    //  const double predicted = line.m_yield*m_rel_act_norms[iso_index(line.m_isotope)];
    //  const double nsigma_off = (predicted - peak.m_counts) / peak.m_counts_uncert;
    //  cout << "\tPeakAt: " << peak.m_energy << " keV (" << line.m_isotope << "), gives predicted counts: "
    //  << predicted << ", verses actual: " << peak.m_counts << " (" << nsigma_off << " sigma diff)\n";
   // }
    
    if( num_isotopes > m_peak_infos.size() )
      throw std::runtime_error( "ManualGenericRelActFunctor: you must have at least as many peaks as"
                               " parameters you are fitting for." );
  }//ManualGenericRelActFunctor constructor
  
  size_t number_residuals() const
  {
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
    return m_peak_infos.size() + 1;
#else
    return m_peak_infos.size();
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
  
  
  double relative_activity( const std::string &iso, const vector<double> &x ) const
  {
    const size_t index = iso_index( iso );
    assert( index < x.size() );
    return m_rel_act_norms[index] * x[index];
  }
  
  
  void eval( const std::vector<double> &x, double *residuals ) const
  {
    m_ncalls += 1;
    
    assert( residuals );
    assert( static_cast<int>(x.size()) == m_isotopes.size() );
    
    assert( m_eqn_order >= 0 );
    const size_t num_eqn_pars = static_cast<size_t>( m_eqn_order + 1 );
    
    const size_t num_isotopes = m_isotopes.size();
    assert( num_isotopes >= 1 );
    
    const double * const pars = x.data();
    vector<double> rel_activities( num_isotopes );
    for( size_t i = 0; i < num_isotopes; ++i )
      rel_activities[i] = this->relative_activity(m_isotopes[i], x);
    
    vector<double> eqn_coefficients;
    vector<vector<double>> eqn_cov;
    RelActCalcManual::fit_rel_eff_eqn_lls( m_eqn_form, m_eqn_order,
                                m_isotopes,
                                rel_activities,
                                m_peak_infos,
                                eqn_coefficients, &eqn_cov );
    
    assert( eqn_coefficients.size() == (m_eqn_order + 1) );
    assert( eqn_cov.size() == (m_eqn_order + 1) );
    
    for( size_t index = 0; index < m_peak_infos.size(); ++index )
    {
      const RelActCalcManual::GenericPeakInfo &peak = m_peak_infos[index];
      
      const double curve_val = eval_eqn( peak.m_energy, m_eqn_form, eqn_coefficients );
      
      double rel_src_counts = 0.0;
      for( const RelActCalcManual::GenericLineInfo &line : peak.m_source_gammas )
      {
        const double rel_activity = relative_activity( line.m_isotope, x );
        rel_src_counts += rel_activity * line.m_yield;
      }//for( const GenericLineInfo &line : peak.m_source_gammas )
      
      
#if( USE_RESIDUAL_TO_BREAK_DEGENERACY )
      if( index == 0 )
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
        if( ((fabs(rel_src_counts) < 1.0E-8) && (rel_src_counts < (1.0E-6*peak.m_counts)))
           || (fabs(rel_src_counts) < std::numeric_limits<float>::epsilon()) )
        {
          rel_src_counts = 1.0E-6 * peak.m_counts;
        }
        
        residuals[index] = (peak.m_counts / rel_src_counts) - curve_val;
      }else if( peak.m_base_rel_eff_uncert == 0.0 )
      {
        // We are not using a m_base_rel_eff_uncert value
        const double pred_counts = curve_val * rel_src_counts;
        residuals[index] = (peak.m_counts - pred_counts) / peak.m_counts_uncert;
      }else
      {
        // We are using a m_base_rel_eff_uncert value
        assert( peak.m_base_rel_eff_uncert <= 1.0 );
        
        const double pred_counts = curve_val * rel_src_counts;
        const double pred_count_uncert = sqrt( pow(peak.m_counts_uncert,2.0)
                                              + pow(rel_src_counts*peak.m_base_rel_eff_uncert,2.0) );
        residuals[index] = (peak.m_counts - pred_counts) / pred_count_uncert;
      }
      
      //cout << "Energy: " << peak.m_energy << " = " << peak.m_source_gammas[0].m_isotope
      //     << " - act=" << relative_activity( peak.m_source_gammas[0].m_isotope, x ) << " :\n";
      //cout << "\trel_efficiency=" << (peak.m_counts / rel_src_counts) << ", curve_val=" << curve_val << endl;
      //cout << "\tpeak.m_counts=" << peak.m_counts << ", rel_src_counts=" << rel_src_counts
      //     << ", (curve_val * rel_src_counts)=" << (curve_val * rel_src_counts) << endl;
      //cout << "\t(rel_efficiency - curve_val):" << ((peak.m_counts / rel_src_counts) - curve_val) << endl;
      //cout << "\t(peak.m_counts - (curve_val * rel_src_counts))/peak.m_counts:" << (peak.m_counts - (curve_val * rel_src_counts))/peak.m_counts << endl;
      //cout << "\t(peak.m_counts - (curve_val * rel_src_counts)) / peak.m_counts_uncert:" << (peak.m_counts - (curve_val * rel_src_counts)) / peak.m_counts_uncert << endl;
      //cout << endl;
    }//for( loop over energies to evaluate at )
    
    //cout << endl << endl << endl;
  }//eval(...)
  
  
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
    for( size_t i = 0; i < m_peak_infos.size(); ++i )
      chi2 += residuals[i]*residuals[i];
    
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
      const vector<double> pars( parameters[0], parameters[0] + m_isotopes.size() );
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

int run_test()
{
  vector<GenericPeakInfo> peak_infos;
  const RelActCalc::RelEffEqnForm eqn_form = RelActCalc::RelEffEqnForm::LnX;
  const size_t eqn_order = 3;  //Three energy dependent terms
  
  GenericPeakInfo info;
  ...
  
  const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
  
  try
  {
    RelEffSolution solution = solve_relative_efficiency( peak_infos, eqn_form, eqn_order );
    
    switch( solution.m_status )
    {
      case ManualSolutionStatus::NotInitialized:
        cout << "Status: NotInitialized\n";
        break;
      case ManualSolutionStatus::ErrorInitializing:
        cout << "Status: ErrorInitializing\n";
        break;
      case ManualSolutionStatus::ErrorFindingSolution:
        cout << "Status: ErrorFindingSolution\n";
        break;
      case ManualSolutionStatus::ErrorGettingSolution:
        cout << "Status: ErrorGettingSolution\n";
        break;
      case ManualSolutionStatus::Success:
        cout << "Status: Success\n";
        break;
    }//switch( solution.m_status )
    
    assert( solution.m_rel_eff_eqn_coefficients.size() == (eqn_order + 1) );
    assert( solution.m_rel_eff_eqn_covariance.size() == (eqn_order + 1) );
    
    cout << "Eqn coefficients: ";
    for( size_t i = 0; i < solution.m_rel_eff_eqn_coefficients.size(); ++i )
    cout << (!i ? "" : ", ") << solution.m_rel_eff_eqn_coefficients[i];
    cout << endl;
    cout << "Eqn coefficient covariance:\n";
    for( size_t i = 0; i < solution.m_rel_eff_eqn_covariance.size(); ++i )
    {
      for( size_t j = 0; j < solution.m_rel_eff_eqn_covariance[i].size(); ++j )
      cout << std::setw(14) << solution.m_rel_eff_eqn_covariance[i][j];
      cout << endl;
    }
    cout << endl;
    
    cout << "Relative activities:" << endl;
    for( const IsotopeRelativeActivity &i : solution.m_rel_activities )
      cout << "\t" << i.m_isotope << ": " << i.m_rel_activity << " +- " << i.m_rel_activity_uncert << endl;
    
    
    const SandiaDecay::SandiaDecayDataBase *db = DecayDataBaseServer::database();
    
    double total_rel_mass = 0.0;
    for( const IsotopeRelativeActivity &i : solution.m_rel_activities )
    {
      const SandiaDecay::Nuclide *nuclide = db->nuclide( i.m_isotope );
      assert( nuclide );
      if( nuclide )
        total_rel_mass += i.m_rel_activity / nuclide->activityPerGram();
    }
    
    cout << "Relative masses:" << endl;
    for( const IsotopeRelativeActivity &i : solution.m_rel_activities )
    {
      const SandiaDecay::Nuclide *nuclide = db->nuclide( i.m_isotope );
      
      if( nuclide )
      {
        const double rel_mass = i.m_rel_activity / nuclide->activityPerGram();
        cout << "\t" << i.m_isotope << ": " << 100.0*rel_mass/total_rel_mass << endl;
      }
    }
      
    
    cout << "Chi2: " << solution.m_chi2 << endl;
    cout << "Num Fcnt Evals: " << solution.m_num_function_eval_total << "\n\n" << endl;
    
    
    for( const auto &peak: peak_infos )
    {
      double function_val = RelActCalc::eval_eqn( peak.m_energy, eqn_form, solution.m_rel_eff_eqn_coefficients );
      
      cout << "For energy " << peak.m_energy << " (";
      for( size_t i = 0; i < peak.m_source_gammas.size(); ++i )
        cout << (i ? ", " : "") << peak.m_source_gammas[i].m_isotope;
      cout << ") function value is " << function_val << ", and:" << endl;
      
      //peak.m_counts
      for( const GenericLineInfo &line : peak.m_source_gammas )
      {
        const string &iso = line.m_isotope;
        const double yield = line.m_yield;
        
        double rel_act = 0.0;
        for( const IsotopeRelativeActivity &r : solution.m_rel_activities )
        {
          if( r.m_isotope == iso )
          {
            rel_act = r.m_rel_activity;
            break;
          }
        }
        
        const double contrib_counts = rel_act * yield * function_val;
        
        cout << "\t" << iso << ": fit " << contrib_counts << " counts and observed "
        << peak.m_counts << "+-" << peak.m_counts_uncert << " (off by "
        << ((contrib_counts - peak.m_counts) / peak.m_counts_uncert)
        << " sigma)" << endl;
      }//for( const RelEff::GammaLineInfo &line : peak.m_source_gammas )
    }//for( const auto &peak: peak_infos )
  }catch( std::exception &e )
  {
    cerr << "Caught exception: " << e.what() << endl;
  }
  
  return EXIT_SUCCESS;
}//int run_test()


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
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  const vector<double> &data_uncertainties = data_uncertainties_orig;
  
  assert( !data_values.empty() );
  assert( energies.size() == data_values.size() );
  assert( energies.size() == data_uncertainties.size() );
  
  if( data_values.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no data points provided." );
  
  const int poly_terms = static_cast<int>(order) + 1;
  const int num_peaks = data_values.size();
  
  
  Eigen::MatrixXd A( num_peaks, poly_terms );
  Eigen::VectorXd b( num_peaks );
  
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const double energy = energies[row];
    const double measured_rel_eff = data_values[row];
    double uncertainty = data_uncertainties[row];
    
    // Letting negative relative efficiencies through doesnt feel right, but I guess we'll do
    //  it to not mess up the L-M fitting of relative activities....
    //if( measured_rel_eff <= 0.0 )
    //  throw runtime_error( "fit_rel_eff_eqn_lls: Measured relative efficiency for energy "
    //                      + to_string(energy) + " is invalid ("
    //                      + to_string(measured_rel_eff) + ")" );
    
    //  But we'll put our foot down for negative or zero uncertainties.
    if( uncertainty <= 0.0 )
      throw runtime_error( "fit_rel_eff_eqn_lls: Uncertainty for energy " + to_string(energy)
                          + " is invalid (" + to_string(uncertainty) + ")" );
    
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
        
        b(row) = std::log(measured_rel_eff) / uncertainty;
        
        break;
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
          A(row,col) = std::pow(std::log(energy), static_cast<double>(col)) / uncertainty;
          break;
          
        case RelActCalc::RelEffEqnForm::LnY:
          //ln(y) = a + b*x + c/x + d/x^2 + e/x^3 + ...
          if( col == 0 )
            A(row,col) = 1.0 / uncertainty;
          else if( col == 1 )
            A(row,col) = energy / uncertainty;
          else
            A(row,col) = std::pow(energy, 1.0 - col) / uncertainty;
          break;
          
        case RelActCalc::RelEffEqnForm::FramEmpirical:
          //ln(y) = a + b/x^2 + c*(lnx) + d*(lnx)^2 + e*(lnx)^3
          if( col == 0 )
            A(row,col) = 1.0 / uncertainty;
          else if( col == 1 )
            A(row,col) = (1.0 / (energy*energy)) / uncertainty;
          else
            A(row,col) = std::pow(std::log(energy), col - 1.0) / uncertainty;
          break;
      }//switch( fcn_form )
    }//for( int col = 0; col < poly_terms; ++col )
  }//for( int col = 0; col < poly_terms; ++col )
  
  // TODO: determine if HouseholderQr or BDC SVD is better/more-stable/faster/whatever
  //const Eigen::VectorXd solution = A.colPivHouseholderQr().solve(b);
  
  const Eigen::BDCSVD<Eigen::MatrixXd> bdc = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::VectorXd solution = bdc.solve(b);
  
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
    const Eigen::MatrixXd A_transpose = A.transpose();
    const Eigen::MatrixXd alpha = Eigen::Product<Eigen::MatrixXd,Eigen::MatrixXd>( A_transpose, A ); //A_transpose * A;
    const Eigen::MatrixXd C = alpha.inverse();
    
    assert( C.rows() == solution.size() );
    assert( C.cols() == solution.size() );
    
    covariance->resize( solution.size() );
    
    for( size_t i = 0; i <= order; ++i )
    {
      vector<double> &row = (*covariance)[i];
      row.resize( solution.size() );
      for( size_t j = 0; j <= order; ++j )
        row[j] = C(i,j);
    }//for( loop over coefficients index )
  }//if( covariance )
}//fit_rel_eff_eqn_lls(...)



void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                         const size_t order,
                         const std::vector<std::string> &isotopes,
                         const std::vector<double> &rel_acts,
                         const std::vector<GenericPeakInfo> &peak_infos,
                         std::vector<double> &fit_pars,
                         std::vector<std::vector<double>> *covariance )
{
  //  We want to solve Ax = b, where
  //    Elements of A are the
  //    x is the coefficients we are solving for
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  
  assert( !isotopes.empty() );
  if( isotopes.empty() )
    throw runtime_error( "fit_rel_eff_eqn_lls: no isotopes specified." );
  
  const int poly_terms = static_cast<int>(order) + 1;
  const int num_peaks = peak_infos.size();
  
  vector<double> energies(num_peaks,0.0), meas_rel_eff(num_peaks,0.0), meas_rel_eff_uncert(num_peaks,0.0);
  
  bool unweighted_fit = false;
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const GenericPeakInfo &peak = peak_infos[row];
    const double energy = peak.m_energy;
    const double counts = peak.m_counts;
    double counts_uncert = peak.m_counts_uncert > 0.0 ? peak.m_counts_uncert : sqrt(peak.m_counts);
    
    double raw_rel_counts = 0.0;
    
    for( const GenericLineInfo &line : peak.m_source_gammas )
    {
      const auto iso_pos = std::lower_bound( std::begin(isotopes), std::end(isotopes), line.m_isotope );
      assert( iso_pos != std::end(isotopes) );
      
      if( iso_pos == std::end(isotopes) )
        throw std::logic_error( "fit_rel_eff_eqn_lls: missing nuclide" );
      
      const size_t iso_index = static_cast<size_t>( iso_pos - std::begin(isotopes) );
      const double rel_act_value = rel_acts[iso_index];
      
      raw_rel_counts += line.m_yield * rel_act_value;
    }//for( const GenericLineInfo &line : peak.m_source_gammas )
    
    
    double measured_rel_eff = counts / raw_rel_counts;
    double measured_rel_eff_uncert = counts_uncert / raw_rel_counts;
    
    // We will clamp rel eff to zero or above ... this is a workaround since Eigens LM doesnt
    //  seem to allow constraining parameter ranges.
    if( (measured_rel_eff <= std::numeric_limits<float>::epsilon())
       || std::isinf(measured_rel_eff)
       || std::isnan(measured_rel_eff) )
    {
      measured_rel_eff = 0.0;
      if( peak.m_base_rel_eff_uncert > std::numeric_limits<float>::epsilon() )
        measured_rel_eff_uncert = 0.0;
      else
        measured_rel_eff_uncert = 1.0;
    }
    
    if( peak.m_base_rel_eff_uncert == -1.0 )
    {
      if( row && !unweighted_fit )
        throw runtime_error( "fit_rel_eff_eqn_lls: for unweighted fit, all peaks must specify m_base_rel_eff_uncert == -1" );
      unweighted_fit = true;
      measured_rel_eff_uncert = 1.0;
    }else
    {
      if( unweighted_fit )
        throw runtime_error( "fit_rel_eff_eqn_lls: for unweighted fit, all peaks must specify m_base_rel_eff_uncert == -1" );
      
      if( (peak.m_base_rel_eff_uncert < 0.0) || (peak.m_base_rel_eff_uncert > 1.0) )
        throw runtime_error( "fit_rel_eff_eqn_lls: m_base_rel_eff_uncert must be in range [0,1]" );
      
      if( peak.m_base_rel_eff_uncert > 0.0 )
      {
        measured_rel_eff_uncert = sqrt( counts_uncert*counts_uncert
                                        + peak.m_base_rel_eff_uncert*peak.m_base_rel_eff_uncert );
      }
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
  
  const double sum_re = std::accumulate( begin(meas_rel_eff), end(meas_rel_eff), 1.0 );
  const double average_re = sum_re / meas_rel_eff.size();
  //const double first_re = meas_rel_eff[0];
  for( double &re : meas_rel_eff )
  {
    re /= average_re;
    //re -= average_re;        // This could go negative
    //re -= (first_re - 1.0);  // Seems to
    //re -= (meas_rel_eff[meas_rel_eff.size()-1] - 1.0);
  }
#endif
  
  fit_rel_eff_eqn_lls( fcn_form, order,
                      energies, meas_rel_eff, meas_rel_eff_uncert,
                      fit_pars, covariance );
}//fit_rel_eff_eqn_lls(...)



void fit_rel_eff_eqn_lls( const RelActCalc::RelEffEqnForm fcn_form,
                         const size_t order,
                         const std::vector<SandiaDecayNucInfo> &nuclides,
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
    GammaInteractionCalc::ShieldingSourceChi2Fcn::cluster_peak_activities( energy_gammas_map,
                                            energy_widths, mixture, n.rel_activity, n.age,
                                            photopeakClusterSigma, energyToCluster, nullptr );
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
    peak.m_energy = energy_gammas[peak_index].first;
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
  
  /*
   // Implementation so we dont have to call into other form of fit_rel_eff_eqn_lls - left commented
   //  out incase changing the form of information somehow ends up actually costing much time.
  const int poly_terms = static_cast<int>(order) + 1;
  const int num_peaks = static_cast<int>(energy_gammas.size());
  
  vector<double> energies(num_peaks,0.0), meas_rel_eff(num_peaks,0.0), meas_rel_eff_uncert(num_peaks,0.0);

  
  bool unweighted_fit = false;
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const double energy = energy_gammas[row].first;
    const double counts = energy_obs_counts[row].second;
    double counts_uncert = energy_obs_counts_uncert[row].second > 0.0 ? energy_obs_counts_uncert[row].second : sqrt(counts);
    
    const double raw_rel_counts = energy_gammas[row].second;
    
    double measured_rel_eff = counts / raw_rel_counts;
    double measured_rel_eff_uncert = counts_uncert / raw_rel_counts;
    
    // We will clamp rel eff to zero or above ... this is a workaround since Eigens LM doesnt
    //  seem to allow constraining parameter ranges.
    if( (measured_rel_eff <= std::numeric_limits<float>::epsilon())
       || std::isinf(measured_rel_eff)
       || std::isnan(measured_rel_eff) )
    {
      measured_rel_eff = 0.0;
      if( base_rel_eff_uncert > std::numeric_limits<float>::epsilon() )
        measured_rel_eff_uncert = 0.0;
      else
        measured_rel_eff_uncert = 1.0;
    }
    
    if( base_rel_eff_uncert == -1.0 )
    {
      unweighted_fit = true;
      measured_rel_eff_uncert = 1.0;
    }else
    {
      if( unweighted_fit )
        throw runtime_error( "fit_rel_eff_eqn_lls: for unweighted fit, all peaks must specify m_base_rel_eff_uncert == -1" );
      
      if( (base_rel_eff_uncert < 0.0) || (base_rel_eff_uncert > 1.0) )
        throw runtime_error( "fit_rel_eff_eqn_lls: m_base_rel_eff_uncert must be in range [0,1]" );
      
      if( base_rel_eff_uncert > 0.0 )
      {
        measured_rel_eff_uncert = sqrt( counts_uncert*counts_uncert + base_rel_eff_uncert*base_rel_eff_uncert );
        //measured_rel_eff_uncert = sqrt( counts_uncert*counts_uncert
        //                                + pow(peak.m_base_rel_eff_uncert*measured_rel_eff, 2.0) );
      }
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
   for( size_t i = 0; i < meas_rel_eff.size(); ++i )
   {
   //  For one test problem doesnt seem to effect the uncertainty on the ratio of activities much
   //  which way of the below we do.  But it doe effect each individual uncertainty a ton (e.g. swaps percentage uncertanty between isotopes)
   meas_rel_eff[i] = meas_rel_eff[i] - (meas_rel_eff[0] - 1.0);
   //meas_rel_eff[i] = meas_rel_eff[i] - (meas_rel_eff[meas_rel_eff.size()-1] - 1.0);
   }
   #endif
   
  fit_rel_eff_eqn_lls( fcn_form, order,
                      energies, meas_rel_eff, meas_rel_eff_uncert,
                      fit_pars, covariance );
   */
}//fit_rel_eff_eqn_lls(...)

void fit_act_to_rel_eff( const RelActCalc::RelEffEqnForm eqn_form,
                        const std::vector<double> &eqn_pars,
                        const std::vector<std::string> &isotopes,
                        const std::vector<GenericPeakInfo> &peak_infos,
                        std::vector<double> &fit_rel_acts,
                        std::vector<double> &fit_uncerts )
{
  // We want to fit the relative activities of the isotopes, given a relative efficiency curve
  
  // We'll start off with some basic sanity checks of the input.
  if( eqn_pars.empty() || eqn_pars.size() > 10 )
    throw runtime_error( "fit_act_to_rel_eff: invalid equation passed in." );
  
  if( peak_infos.size() < isotopes.size() )
    throw runtime_error( "fit_act_to_rel_eff: less peaks than isotopes." );
  
  std::vector<bool> used_isotopes( isotopes.size(), false );
  for( const GenericPeakInfo &info : peak_infos )
  {
    // We could blissfully ignore peaks with no source gammas (and assign them an activity of zero),
    //  but we'll be strict to maybe prevent some input mistakes, or something.
    if( info.m_source_gammas.empty() )
      throw runtime_error( "fit_act_to_rel_eff: peaks at "
                          + to_string(info.m_energy) + " keV has no source gammas." );
    
    for( const GenericLineInfo &line : info.m_source_gammas )
    {
      if( line.m_yield <= 0.0 || std::isinf(line.m_yield) || std::isnan(line.m_yield) )
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
  }//for( const RelEff::PeakInfo &info : peak_infos )
  
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
  
  
  const int num_isotopes = static_cast<int>(isotopes.size());
  const int num_peaks = peak_infos.size();
  
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero( num_peaks, num_isotopes );
  Eigen::VectorXd b( num_peaks );
  
  
  for( size_t row = 0; row < num_peaks; ++row )
  {
    const GenericPeakInfo &peak = peak_infos[row];
    const double energy = peak.m_energy;
    const double counts = peak.m_counts;
    const double counts_uncert = peak.m_counts_uncert;
    
    const double rel_eff = eval_eqn( energy, eqn_form, eqn_pars );
    
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
  
  const Eigen::BDCSVD<Eigen::MatrixXd> bdc = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
  const Eigen::VectorXd solution = bdc.solve(b);
  
  assert( solution.size() == num_isotopes );
  
  // TODO: I'm sure Eigen::BDCSVD has the uncertainty matrix in it somewhere already computed, but
  //       for the moment (See pg 796 in Numerical Recipes for hint - probably has something to do
  //       with bdc.singularValues(), bdc.matrixV(), or bdc.matrixU()) we'll just be dumb and do
  //       extra (unstable?) work.
  const Eigen::MatrixXd A_transpose = A.transpose();
  const Eigen::MatrixXd alpha = Eigen::Product<Eigen::MatrixXd,Eigen::MatrixXd>( A_transpose, A ); //A_transpose * A;
  const Eigen::MatrixXd C = alpha.inverse();
  
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
  return eval_eqn( energy, m_rel_eff_eqn_form, m_rel_eff_eqn_coefficients );
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


double RelEffSolution::ratio( const size_t iso1, const size_t iso2 ) const
{
  assert( iso1 < m_rel_activities.size() );
  assert( iso2 < m_rel_activities.size() );
  return m_rel_activities[iso1].m_rel_activity / m_rel_activities[iso2].m_rel_activity;
}

double RelEffSolution::ratio( const std::string &nuc1, const std::string &nuc2 ) const
{
  return ratio( nuclide_index(nuc1), nuclide_index(nuc2) );
}


double RelEffSolution::ratio_uncert( const size_t iso1, const size_t iso2 ) const
{
  assert( iso1 != iso2 );
  assert( iso1 < m_rel_act_covariance.size() );
  assert( iso2 < m_rel_act_covariance.size() );
  assert( m_rel_act_covariance.size() == m_activity_norms.size() );
  
  const double norm_1 = m_activity_norms[iso1];
  const double norm_2 = m_activity_norms[iso2];
  
  // Divide the relative activities by their norms, to go back to the actual values fit for.
  const double fit_act_1 = m_rel_activities[iso1].m_rel_activity / norm_1;
  const double fit_act_2 = m_rel_activities[iso2].m_rel_activity / norm_2;
  
  const double cov_1_1 = m_rel_act_covariance[iso1][iso1];
  const double cov_1_2 = m_rel_act_covariance[iso1][iso2];
  const double cov_2_2 = m_rel_act_covariance[iso2][iso2];
  
  
  const double ratio = m_rel_activities[iso1].m_rel_activity
  / m_rel_activities[iso2].m_rel_activity;
  
  // TODO: I think this is the right way to compute ratio, taking into account correlations, given
  //       we actually fit for the values that multiplied m_activity_norms[],... need to double check
  const double uncert = fabs(ratio)
  * sqrt( (cov_1_1 / fit_act_1 / fit_act_1)
         + (cov_2_2 / fit_act_2 / fit_act_2)
         - (2.0 * cov_1_2 / fit_act_1 / fit_act_2) );
  
  cout
  << "Ratio " << m_rel_activities[iso1].m_isotope << "/" << m_rel_activities[iso2].m_isotope
  << " = " << ratio << " +- " << uncert << " (would be +- "
  << fabs(ratio) * sqrt( m_rel_activities[iso1].m_rel_activity_uncert*m_rel_activities[iso1].m_rel_activity_uncert
                        + m_rel_activities[iso2].m_rel_activity_uncert*m_rel_activities[iso2].m_rel_activity_uncert )
  << " w/ No correlation)" << endl;
  
  return uncert;
}

double RelEffSolution::ratio_uncert( const std::string &iso1, const std::string &iso2 ) const
{
  return ratio_uncert( nuclide_index(iso1), nuclide_index(iso2) );
}



RelEffSolution solve_relative_efficiency( const std::vector<GenericPeakInfo> &peak_infos,
                                          const RelActCalc::RelEffEqnForm eqn_form,
                                          const size_t eqn_order )
{
  const auto start_time = std::chrono::high_resolution_clock::now();
  
  RelEffSolution solution;
  
  DoWorkOnDestruct setFinalTime( [&solution,start_time](){
    const auto end_time = std::chrono::high_resolution_clock::now();
    solution.m_num_microseconds_eval = std::chrono::duration<double, std::micro>(end_time - start_time).count();
  });
  
  solution.m_rel_eff_eqn_form = eqn_form;
  solution.m_rel_eff_eqn_order = eqn_order;
  solution.m_status = ManualSolutionStatus::NotInitialized;
  solution.m_input_peak = peak_infos;
  
  ManualGenericRelActFunctor *cost_functor = nullptr;
  try
  {
    cost_functor = new ManualGenericRelActFunctor( eqn_form, eqn_order, peak_infos );
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
  
  const size_t num_peaks = cost_functor->m_peak_infos.size();
  const size_t num_nuclides = cost_functor->m_isotopes.size();
  
  solution.m_activity_norms = cost_functor->m_rel_act_norms;
  
  auto cost_function = new ceres::DynamicNumericDiffCostFunction<ManualGenericRelActFunctor>( cost_functor );
  // We have one more residual than the number of peaks; the last residual is what clamps the
  //  relative efficiency curve to 1.0 at the lowest energy (to prevent degeneracy between rel eff
  //  curve, and rel activities).
  cost_function->SetNumResiduals( cost_functor->number_residuals() );
  cost_function->AddParameterBlock( num_nuclides );
  
  // Relative activities multiples start out as 1.0 because ManualGenericRelActFunctor constructor
  //   estimates the activities for a flat rel eff = 1.0; see
  //   #ManualGenericRelActFunctor::m_rel_act_norms.
  vector<double> parameters( num_nuclides, 1.0 );
  double *pars = &parameters[0];
  
  vector<double *> parameter_blocks( 1, pars );
  
  ceres::Problem problem;
  
  // TODO: investigate using a LossFunction - probably really need it
  ceres::LossFunction *lossfcn = nullptr;
  problem.AddResidualBlock( cost_function, lossfcn, parameter_blocks );
  
  // Set a lower bound on relative activities to be 0
  for( size_t i = 0; i < num_nuclides; ++i )
    problem.SetParameterLowerBound( pars, i, 0.0 );
  
  
  // Okay - we've set our problem up
  ceres::Solver::Options ceres_options;
  ceres_options.linear_solver_type = ceres::DENSE_QR;
  ceres_options.minimizer_progress_to_stdout = true; //true;
  
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
  
  
  ceres::Solver::Summary summary;
  
  try
  {
    ceres::Solve(ceres_options, &problem, &summary);
    
    solution.m_num_function_eval_solution = cost_functor->m_ncalls;
    
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
  
  
  //std::cout << summary.BriefReport() << "\n";
  std::cout << summary.FullReport() << "\n";
  cout << "Took " << solution.m_num_function_eval_solution << " calls to solve." << endl;
  
  cout << "\nFit rel_act parameters: {";
  for( size_t i = 0; i < parameters.size(); ++i )
  cout << (i ? ", " : "") << parameters[i];
  cout << "}\n\n";
  
  //cout << "\n\n\nHackign in values of parameters!!!\n\n" << endl;
  //parameters[0] = 0.528839;
  //parameters[1] = 1.2962;
  //parameters[2] = 0.538641;
  
  try
  {
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::CovarianceAlgorithmType::SPARSE_QR; //
    cov_options.num_threads = ceres_options.num_threads;
    
    vector<double> uncertainties( num_nuclides, 0.0 ), uncerts_squared( num_nuclides, 0.0 );
    
    ceres::Covariance covariance(cov_options);
    vector<pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back( {pars, pars} );
    
    solution.m_rel_act_covariance.clear();
    if( !covariance.Compute(covariance_blocks, &problem) )
    {
      cerr << "Failed to compute final covariances!" << endl;
      solution.m_warnings.push_back( "Failed to compute final covariances." );
    }else
    {
      // row-major order: the elements of the first row are consecutively given, followed by second
      //                  row contents, etc
      vector<double> row_major_covariance( num_nuclides * num_nuclides );
      covariance.GetCovarianceBlock( pars, pars, row_major_covariance.data() );
      
      solution.m_rel_act_covariance.resize( num_nuclides, vector<double>(num_nuclides,0.0) );
      for( size_t row = 0; row < num_nuclides; ++row )
      {
        for( size_t col = 0; col < num_nuclides; ++col )
        solution.m_rel_act_covariance[row][col] = row_major_covariance[row*num_nuclides + col];
      }//for( size_t row = 0; row < num_nuclides; ++row )
    }//if( we failed to get covariance ) / else
    
    
    vector<double> rel_activities( num_nuclides );
    
    for( size_t i = 0; i < num_nuclides; ++i )
    rel_activities[i] = cost_functor->relative_activity( cost_functor->m_isotopes[i], parameters );
    
    fit_rel_eff_eqn_lls( eqn_form, eqn_order, cost_functor->m_isotopes, rel_activities, peak_infos,
                        solution.m_rel_eff_eqn_coefficients, &(solution.m_rel_eff_eqn_covariance) );
    
    assert( solution.m_rel_eff_eqn_coefficients.size() == (eqn_order + 1) );
    
    for( size_t i = 0; i < cost_functor->m_isotopes.size(); ++i )
    {
      const double rel_act_norm = cost_functor->m_rel_act_norms[i];
      
      IsotopeRelativeActivity rel_act;
      rel_act.m_isotope = cost_functor->m_isotopes[i];
      
      rel_act.m_rel_activity = rel_act_norm * parameters[i];
      if( solution.m_rel_act_covariance.empty() )
      {
        rel_act.m_rel_activity_uncert = -1.0;
      }else
      {
        assert( i < solution.m_rel_act_covariance.size() );
        assert( i < solution.m_rel_act_covariance[i].size() );
        rel_act.m_rel_activity_uncert = rel_act_norm * std::sqrt( solution.m_rel_act_covariance[i][i] );
      }
      
      solution.m_rel_activities.push_back( std::move(rel_act) );
    }//for( loop over relative activities )
    
    
    // We'll manually compute the Chi2 here, so we only take into account statistical uncertainties,
    //  (i.e., we ignore #GenericPeakInfo::m_base_rel_eff_uncert)
    solution.m_chi2 = 0.0;
    
    // TODO: The DOF is probably off by one - need to think on this and come back to it
    assert( cost_functor->m_peak_infos.size() >= ((eqn_order+1) + (cost_functor->m_isotopes.size() - 1)) );
    solution.m_dof = cost_functor->m_peak_infos.size() - (eqn_order+1) - (cost_functor->m_isotopes.size() - 1);
    
    for( const GenericPeakInfo &peak : cost_functor->m_peak_infos )
    {
      const vector<double> &rel_eff_coefs = solution.m_rel_eff_eqn_coefficients;
      const double curve_val = RelActCalc::eval_eqn( peak.m_energy, eqn_form, rel_eff_coefs );
      
      double expected_src_counts = 0.0;
      for( const GenericLineInfo &line : peak.m_source_gammas )
      {
        const double rel_activity = cost_functor->relative_activity( line.m_isotope, parameters );
        expected_src_counts += rel_activity * line.m_yield;
      }//for( const RelEff::GammaLineInfo &line : peak.m_source_gammas )
      
      const double expected_counts = expected_src_counts * curve_val;
      solution.m_chi2 += std::pow( (expected_counts - peak.m_counts) / peak.m_counts_uncert, 2.0 );
    }//for( loop over energies to evaluate at )
    
    
    solution.m_status = ManualSolutionStatus::Success;
    solution.m_num_function_eval_total = cost_functor->m_ncalls;
  }catch( std::exception &e )
  {
    solution.m_status = ManualSolutionStatus::ErrorGettingSolution;
    solution.m_error_message = e.what();
    cerr << "RelActCalcManual::solve_relative_efficiency: Failed to get solution after solving: "
         << e.what() << endl;
    
    return solution;
  }//try / catch to get the solution
  
  
  return solution;
}//solve_relative_efficiency(...)

}//namespace RelActCalcManual
