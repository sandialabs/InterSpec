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

#include <vector>

#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MinosError.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/MakeDrfFit.h"


using namespace std;



namespace
{
  class DetectorResolutionFitness
  : public ROOT::Minuit2::FCNBase
  {
  protected:
    std::deque<std::shared_ptr<const PeakDef>> m_peaks;
    DetectorPeakResponse::ResolutionFnctForm m_form;
    
  public:
    DetectorResolutionFitness( const std::deque<std::shared_ptr<const PeakDef>> &peaks,
                              DetectorPeakResponse::ResolutionFnctForm form )
    : m_peaks( peaks ),
    m_form( form )
    {
    }
    
    virtual ~DetectorResolutionFitness(){}
    
    
    virtual double DoEval( const std::vector<double> &x ) const
    {
      for( size_t i = 0; i < x.size(); ++i )
        if( IsInf(x[i]) || IsNan(x[i]) )
          return 999999.9;
      
      vector<float> fx( x.size() );
      for( size_t i = 0; i < x.size(); ++i )
        fx[i] = static_cast<float>( x[i] );
      
      double chi2 = 0.0;
      for( const PeakModel::PeakShrdPtr &peak : m_peaks )
      {
        if( !peak )  //probably isnt needed
          continue;
        const float predicted = DetectorPeakResponse::peakResolutionSigma( peak->mean(), m_form, fx );
        if( predicted <= 0.0 || IsNan(predicted) || IsInf(predicted) )
          return 999999.9;
        
        chi2 += MakeDrfFit::peak_width_chi2( predicted, *peak );
      }//for( const EnergySigma &es : m_energies_and_sigmas )
      
      return chi2;
    }//DoEval();
    
    virtual double operator()( const std::vector<double> &x ) const
    {
      return DoEval( x );
    }
    
    DetectorResolutionFitness &operator=( const DetectorResolutionFitness &rhs )
    {
      if( &rhs != this )
      {
        m_peaks = rhs.m_peaks;
        m_form = rhs.m_form;
      }
      return *this;
    }
    
    virtual double Up() const{ return 1.0; }
  };//class DetectorResolutionFitness
  
  
  class DetectorEffFitness
  : public ROOT::Minuit2::FCNBase
  {
  protected:
    std::vector<MakeDrfFit::DetEffDataPoint> m_data;
    int m_order;
    
  public:
    DetectorEffFitness( const std::vector<MakeDrfFit::DetEffDataPoint> &data, const int order )
    : m_data( data ),
      m_order( order )
    {
    }
    
    virtual ~DetectorEffFitness(){}
    
    virtual double DoEval( const std::vector<double> &x ) const
    {
      for( const double val : x )
        if( IsInf(val) || IsNan(val) )
          return 999999.9;
      
      vector<float> fx( x.size() );
      for( size_t i = 0; i < x.size(); ++i )
        fx[i] = static_cast<float>( x[i] );
      
      double chi2 = 0.0;
      for( const MakeDrfFit::DetEffDataPoint &data : m_data )
      {
        const float eqneff = DetectorPeakResponse::expOfLogPowerSeriesEfficiency( data.energy, fx );
        if( eqneff <= 0.0 || IsNan(eqneff) || IsInf(eqneff) )
          return 999999.9;
        
        const double uncert = data.efficiency_uncert <= 0.0 ? 0.05*data.efficiency : data.efficiency_uncert;
        //cout << "data.energy=" << data.energy << ", eqneff=" << eqneff
        //     << ", data.efficiency=" << data.efficiency << ", uncert=" << uncert
        //     << ", data.efficiency_uncert=" << data.efficiency_uncert << endl;
        chi2 += std::pow( (eqneff - data.efficiency)/uncert, 2 );
      }
      
      return chi2;
    }//DoEval();
    
    virtual double operator()( const std::vector<double> &x ) const
    {
      return DoEval( x );
    }
    
    DetectorEffFitness &operator=( const DetectorEffFitness &rhs )
    {
      if( &rhs != this )
      {
        m_data = rhs.m_data;
        m_order = rhs.m_order;
      }
      return *this;
    }
    
    virtual double Up() const{ return 1.0; }
  };//class DetectorResolutionFitness
}//namespace




namespace MakeDrfFit
{
  
double peak_width_chi2( double predicted_sigma, const PeakDef &peak )
{
  if( !peak.gausPeak() )
    return 0.0;  //shouldnt ever get here
  
  const double measured_sigma = peak.sigma();
  //double measured_sigma_uncert = peak.sigmaUncert();
  //if( measured_sigma_uncert <= 0.01 )
  //  measured_sigma_uncert = 1.0;
  const double measured_sigma_uncert = ((peak.sigmaUncert() > 0.0) ? std::max( peak.sigmaUncert(), 0.01*measured_sigma) : 0.05*measured_sigma);
  
  const double chi = (measured_sigma - predicted_sigma)/measured_sigma_uncert;
  return chi*chi;
}//peak_width_chi2(...)



double performResolutionFit( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                           const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                            const bool highres,
                           const int sqrtEqnOrder,
                           std::vector<float> &answer,
                           std::vector<float> &uncerts )
{
  if( !peaks || peaks->empty() )
    throw runtime_error( "MakeDrfFit::performResolutionFit(...): no input peaks" );
  
  bool fit_using_lls = false;
  double a_initial, b_initial, c_initial;
  double lowerA, upperA, lowerB, upperB, lowerC, upperC;
  
  switch( fnctnlForm )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
    {
      if( highres )
      {
        lowerA = 0.75*1.0;
        upperA = 2.0*1.77;
        lowerB = 0.75*0.20028;
        upperB = 2.0*0.27759;
        lowerC = 0.75*0.31;
        upperC = 1.5*0.57223;
      }else
      {
        lowerA = 1.5*-7.0;
        upperA = 1.5*7.44000;
        lowerB = 0.75*2.13000;
        upperB = 1.5*8.50000;
        lowerC = 0.75*0.20000;
        upperC = 1.25*0.70000;
      }
      
      if( answer.size() == 3 )
      { //caller has provided some default values, lets use them
        a_initial = static_cast<double>( answer[0] );
        b_initial = static_cast<double>( answer[1] );
        c_initial = static_cast<double>( answer[2] );
      }else
      {
        if( highres )
        {
          a_initial = 0.5*(1.0 + 1.77);
          b_initial = 0.5*(0.20028 + 0.27759);
          c_initial = 0.5*(0.31 + 0.57223);
        }else
        {
          a_initial = 0.0;//0.5*(-7.0 + 7.44000);
          b_initial = 0.5*(2.13000 + 8.50000);
          c_initial = 0.5*(0.20000 + 0.70000);
        }
      }//if( answer.size() == 3 ) / else
      
      break;
    }//case kGadrasResolutionFcn:
      
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
    {
      if( highres )
      {
        //Based on pretty much nothing
        a_initial = 2.6;
        lowerA = -2.5;
        upperA = 7.5;
        
        b_initial = 1.0;
        lowerB = -5.0;
        upperB = 5.0;
        
        c_initial = 0;
        lowerC = -5.0;
        upperC = 5.0;
      }else
      {
        //Based on zero detectors so far
        a_initial = 100;
        lowerA = -100;
        upperA = 400;
        
        b_initial = 3600.0 / 0.661;
        lowerB = 0.0;
        upperB = 180*180 / 0.661;
        
        c_initial = 0;
        lowerC = -10000.0;
        upperC = 10000.0;
      }//if( highres ) / else
      
      try
      {
        int num_fit_coefficients = 3;
        //if( sqrtEqnOrder > 3 )
        //  num_fit_coefficients = sqrtEqnOrder;  // Right now, just assuming 3 fit parameters always, but thing would work out if we wanted to do different
        double chi2 = MakeDrfFit::fit_sqrt_poly_fwhm_lls( *peaks, num_fit_coefficients, true, answer, uncerts );
        
        assert( answer.size() == static_cast<size_t>(num_fit_coefficients) );
        
        //cout << "MakeDrfFit::fit_sqrt_poly_fwhm_lls got {";
        //for( size_t i = 0; i < answer.size(); ++i )
        //  cout << answer[i] << "+-" << uncerts[i] << ", ";
        //cout << "}.  Chi2=" << chi2 << endl;
        
        fit_using_lls = true;
      }catch( std::exception &e )
      {
        cerr << "MakeDrfFit::fit_sqrt_poly_fwhm_lls threw exception: " << e.what() << endl;
        
      }//try / catch
      
      break;
    }//case DetectorPeakResponse::kSqrtEnergyPlusInverse:
      
      
    case DetectorPeakResponse::kSqrtPolynomial:
    {
      if( sqrtEqnOrder < 1 )
        throw runtime_error( "performResolutionFit: sqrt eqn order should be at least 1" );
        
      if( highres )
      {
        //Based on pretty much nothing
        a_initial = 2.6;
        lowerA = -2.5;
        upperA = 7.5;
        
        b_initial = 1.0;
        lowerB = -5.0;
        upperB = 5.0;
        
        c_initial = 0;
        lowerC = -5.0;
        upperC = 5.0;
      }else
      {
        //Based on zero detectors so far
        a_initial = 100;
        lowerA = -100;
        upperA = 400;
        
        b_initial = 3600.0 / 0.661;
        lowerB = 0.0;
        upperB = 180*180 / 0.661;
        
        c_initial = 0;
        lowerC = -3600.0/(0.661*0.661);
        upperC = 3600.0/(0.661*0.661);
      }//if( highres ) / else
      
      try
      {
        double chi2 = MakeDrfFit::fit_sqrt_poly_fwhm_lls( *peaks, sqrtEqnOrder, false, answer, uncerts );
        
        assert( answer.size() == static_cast<int>(sqrtEqnOrder) );
        
        //cout << "MakeDrfFit::fit_sqrt_poly_fwhm_lls got {";
        //for( size_t i = 0; i < answer.size(); ++i )
        //  cout << answer[i] << "+-" << uncerts[i] << ", ";
        //cout << "}.  Chi2=" << chi2 << endl;
        
        fit_using_lls = true;
      }catch( std::exception &e )
      {
        cerr << "MakeDrfFit::fit_sqrt_poly_fwhm_lls threw exception: " << e.what() << endl;
        
      }//try / catch
      
      break;
    }//case kSqrtPolynomial:
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
      //We should never make it here if the detector FWHM functional form has
      //  not been explicitly set, so throw and error if we got here
      throw runtime_error( "MakeDrfFit::performResolutionFit(...):"
                          " invalid ResolutionFnctForm" );
      break;
  }//switch( fnctnlForm )
  
  
  //For high res we should only have to to this one fit
  //But for lowres detectors we should consider the cases of
  //  A==0, A<0, and A>0 when fitting for kGadrasResolutionFcn
  
  bool fitting_gad_a = false;
  DetectorResolutionFitness fitness( *peaks, fnctnlForm );
  
  ROOT::Minuit2::MnUserParameters inputPrams;
  switch( fnctnlForm )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
    {
      if( peaks->size() == 1 )
      {
        inputPrams.Add( "A", a_initial );
        inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
        inputPrams.Add( "C", c_initial );
      }else if( peaks->size() == 2 )
      {
        inputPrams.Add( "A", a_initial );
        inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
        inputPrams.Add( "C", c_initial, 0.1*(upperC-lowerC), lowerC, upperC );
      }else
      {
        fitting_gad_a = true;
        inputPrams.Add( "A", a_initial, 0.1*(upperA-lowerA), lowerA, upperA );
        inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
        inputPrams.Add( "C", c_initial, 0.1*(upperC-lowerC), lowerC, upperC );
      }//if( all_peaks.size() == 1 )
      break;
    }//case DetectorPeakResponse::kGadrasResolutionFcn:
      
      
    case DetectorPeakResponse::kSqrtEnergyPlusInverse:
    {
      if( fit_using_lls )
      {
        for( int order = 0; order < 3; ++order )
        {
          const string name = string("") + char('A' + order);
          inputPrams.Add( name, answer[order], uncerts[order] );
        }//for( int order = 0; order < sqrtEqnOrder; ++order )
      }else
      {
        if( peaks->size() == 1 )
        {
          inputPrams.Add( "A", 0.0 );
          inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
          inputPrams.Add( "C", 0.0 );
        }if( peaks->size() == 2 )
        {
          inputPrams.Add( "A", a_initial, 0.1*(upperA-lowerA), lowerA, upperA );
          inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
          inputPrams.Add( "C", 0.0 );
        }else if( peaks->size() >= 3 )
        {
          inputPrams.Add( "A", a_initial, 0.1*(upperA-lowerA), lowerA, upperA );
          inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
          inputPrams.Add( "C", c_initial, 0.1*(upperC-lowerC), lowerC, upperC );
        }
      }//if( fit_using_lls )
      
      break;
    }//case DetectorPeakResponse::kSqrtEnergyPlusInverse:
      
      
    case DetectorPeakResponse::kSqrtPolynomial:
    {
      if( fit_using_lls )
      {
        for( int order = 0; order < sqrtEqnOrder; ++order )
        {
          const string name = string("") + char('A' + order);
          inputPrams.Add( name, answer[order], uncerts[order] );
        }//for( int order = 0; order < sqrtEqnOrder; ++order )
      }else
      {
        if( peaks->size() < 3 )
        {
          inputPrams.Add( "A", a_initial, 0.1*(upperA-lowerA), lowerA, upperA );
          inputPrams.Add( "B", b_initial );
        }else
        {
          inputPrams.Add( "A", a_initial, 0.1*(upperA-lowerA), lowerA, upperA );
          inputPrams.Add( "B", b_initial, 0.1*(upperB-lowerB), lowerB, upperB );
          if( peaks->size() > 3 )
            inputPrams.Add( "C", c_initial, 0.1*(upperC-lowerC), lowerC, upperC );
          //inputPrams.Add( "D", 0.5, 0.05, 0.25, 0.75 );
        }
      }//if( fit_using_lls )
      break;
    }//case kSqrtPolynomial:
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
      assert( 0 );
      break;
  }//switch
  
  
  ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( fitness, inputParamState, strategy );
  
  double tolerance = 0.5;
  unsigned int maxFcnCall = 50000;
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
  
  ROOT::Minuit2::MnUserParameters fitParams = fitter.Parameters();
  
  auto scanpar = [&fitParams,&minimum,&fitness,&fitter,maxFcnCall,tolerance,strategy]( const string par ){
    std::vector<double> best_pars = fitter.Params();
    
    const auto index = fitter.Parameters().Index(par);
    auto &mnpar = fitter.Parameter( index );
    
    assert( mnpar.HasLowerLimit() && mnpar.HasUpperLimit() );
    
    double best_a = mnpar.Value();
    const double lower_a = mnpar.LowerLimit();
    const double upper_a = mnpar.UpperLimit();
    
    double best_chi2 = fitness.DoEval( fitParams.Params() );
    
    for( double an = lower_a; an < upper_a; an += 0.1*(upper_a-lower_a) )
    {
      auto testpar = fitParams;
      testpar.SetValue(par, an);
      testpar.Fix(par);
      
      ROOT::Minuit2::MnUserParameterState anInputParam( testpar );
      ROOT::Minuit2::MnMinimize anfitter( fitness, anInputParam, strategy );
      
      ROOT::Minuit2::FunctionMinimum anminimum = anfitter( maxFcnCall, tolerance );
      if( !anminimum.IsValid() )
        anminimum = anfitter( maxFcnCall, tolerance );
      
      const double this_chi2 = fitness.DoEval( anfitter.Params() );
      if( this_chi2 < best_chi2 )
      {
        best_pars = anfitter.Params();
        best_chi2 = this_chi2;
        best_a = an;
      }
    }//for( double an = 1.0; an < 101.0; an += 2.5 )
    
    if( best_chi2 < fitness.DoEval(fitParams.Params()) )
    {
      auto testpar = fitParams;
      testpar.SetValue( par, best_pars[index] );
      ROOT::Minuit2::MnUserParameterState anInputParam( testpar );
      ROOT::Minuit2::MnMinimize anfitter( fitness, anInputParam, strategy );
      
      minimum = anfitter( maxFcnCall, tolerance );
      if( !minimum.IsValid() )
        minimum = anfitter( maxFcnCall, tolerance );
      
      fitParams = anfitter.Parameters();
    }//if( best_chi2 < fitness.DoEval(fitParams.Params()) )
  };//scanpar lambda
  
  if( fitting_gad_a )
  {
    scanpar( "A" );
  }//if( fitting_gad_a )
  
  const double final_chi2 = fitness.DoEval( fitParams.Params() );
  
  //cout << "FWHM final chi2=" << final_chi2 << endl;
  if( fit_using_lls )
  {
    const double pre_chi2 = fitness.DoEval( vector<double>( begin(answer), end(answer) ) );
    //cout << "FWHM LLS chi2=" << pre_chi2 << endl;
    if( pre_chi2 < final_chi2 )
    {
      //cout << "Least Linear chi2 better than from Minuit, using that" << endl;
      return pre_chi2;
    }
  }
  
  if( !minimum.IsValid() )
  {
    stringstream msg;
    msg  << "FWHM response function fit status is not valid"
    << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
    << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
    << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
    << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
    << "\n\tHasValidParameters: " << minimum.HasValidParameters()
    << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
    << endl;
    if( minimum.IsAboveMaxEdm() )
      msg << "\t\tEDM=" << minimum.Edm() << endl;
    cerr << endl << msg.str() << endl;
    throw std::runtime_error( msg.str() );
  }//if( !minimum.IsValid() )
  
  answer.clear();
  uncerts.clear();
  
  //for( size_t i = 0; i < fitParams.Params().size(); ++i )
  //  cout << "\tMinuit fit FWHM Par_" << i << "=" << fitParams.Params()[i] << endl;
  
  for( const double p : fitParams.Params() )
    answer.push_back( static_cast<float>(p) );
  
  for( const double p : fitParams.Errors() )
    uncerts.push_back( p );
  
  return final_chi2;
}//std::vector<float> performResolutionFit(...)


template<class T>
bool matrix_invert( const boost::numeric::ublas::matrix<T>& input,
                     boost::numeric::ublas::matrix<T> &inverse )
{
  using namespace boost::numeric;
  ublas::matrix<T> A( input );
  ublas::permutation_matrix<std::size_t> pm( A.size1() );
  const size_t res = lu_factorize(A, pm);
  if( res != 0 )
    return false;
  inverse.assign( ublas::identity_matrix<T>( A.size1() ) );
  lu_substitute(A, pm, inverse);
  return true;
}//matrix_invert
  
  
double fit_sqrt_poly_fwhm_lls( const std::deque< std::shared_ptr<const PeakDef> > &peaks,
                               const int num_fit_coefficients,
                               const bool include_inv_term,
                               std::vector<float> &coeffs,
                               std::vector<float> &coeff_uncerts )
{
  const size_t nbin = peaks.size();
  
  if( num_fit_coefficients < 1 )
    throw runtime_error( "fit_sqrt_poly_fwhm_lls: num_fit_coefficients must be >= 1" );
  
  if( !nbin )
    throw runtime_error( "fit_sqrt_poly_fwhm_lls: must have at least 1 input peak" );
  
  if( nbin < static_cast<size_t>(num_fit_coefficients) )
    throw runtime_error( "fit_sqrt_poly_fwhm_lls: must have at least as many peaks as coefficients to fit to" );
  
  //log(eff(x)) = A0 + A1*logx + A2*logx^2 + A3*logx^3, where x is energy in MeV
  vector<float> x, widths, widths_uncert;
  x.resize( peaks.size() );
  widths.resize( peaks.size() );
  widths_uncert.resize( peaks.size() );
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    if( peaks[i]->gausPeak() )
    {
      x[i] = peaks[i]->mean() / (include_inv_term ? 1.0f : 1000.0f);
      widths[i] = peaks[i]->fwhm();
      widths_uncert[i] = 2.35482*((peaks[i]->sigmaUncert() > 0.0) ? std::max( peaks[i]->sigmaUncert(), 0.01*widths[i]) : 0.05*widths[i]);
    }
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  //General Linear Least Squares fit
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation is quite inneficient.
  using namespace boost::numeric;
  
  ublas::matrix<double> A( nbin, num_fit_coefficients );
  ublas::vector<double> b( nbin );
  
  
  for( size_t row = 0; row < nbin; ++row )
  {
    const double data_y = widths[row] * widths[row];
    const double data_y_uncert = widths_uncert[row] * widths_uncert[row];
    
    b(row) = data_y / data_y_uncert;
    
    if( include_inv_term )
    {
      // Nominally we will only ever call this function to fit for the FRAM style FWHM  when num_fit_coefficients == 3,
      //  but we'll go ahead and code in possibility to fit for different num_fit_coefficients.  For num_fit_coefficients==1, we'll
      //  fit the energy dependent term, for num_fit_coefficients==2, the constant + energy dependent, and num_fit_coefficients==3
      //  the inverse term as well.  Higher num_fit_coefficients terms we'll do as power series in energy.
      //  Before we return, we'll swap the zeroth and first coefficients so things are as expected.
      assert( num_fit_coefficients == 3 );
      for( int col = 0; col < num_fit_coefficients; ++col )
      {
        if( col == 0 )
          A(row,col) = x[row] / data_y_uncert;
        else if( col == 1 )
          A(row,col) = 1.0 / data_y_uncert;
        else if( col == 2 )
          A(row,col) = (1.0/x[row]) / data_y_uncert;
        else
          A(row,col) = std::pow( x[row], double(col-1)) / data_y_uncert;
      }
    }else
    {
      for( int col = 0; col < num_fit_coefficients; ++col )
        A(row,col) = std::pow( x[row], double(col)) / data_y_uncert;
    }
  }//for( int col = 0; col < num_fit_coefficients; ++col )
  
  const ublas::matrix<double> A_transpose = ublas::trans( A );
  const ublas::matrix<double> alpha = prod( A_transpose, A );
  ublas::matrix<double> C( alpha.size1(), alpha.size2() );
  const bool success = matrix_invert( alpha, C );
  if( !success )
    throw runtime_error( "fit_sqrt_poly_fwhm_lls(...): trouble inverting matrix" );
  
  const ublas::vector<double> beta = prod( A_transpose, b );
  const ublas::vector<double> a = prod( C, beta );
  
  coeffs.resize( num_fit_coefficients );
  coeff_uncerts.resize( num_fit_coefficients );
  for( int coef = 0; coef < num_fit_coefficients; ++coef )
  {
    coeffs[coef] = static_cast<float>( a(coef) );
    coeff_uncerts[coef] = static_cast<float>( C(coef,coef) );
  }//for( int coef = 0; coef < num_fit_coefficients; ++coef )
  
  double chi2 = 0;
  for( size_t bin = 0; bin < nbin; ++bin )
  {
    double y_pred = 0.0;
    if( include_inv_term )
    {
      for( int i = 0; i < num_fit_coefficients; ++i )
      {
        if( i == 0 )
          y_pred += a(i) * x[bin];
        else if( i == 1 )
          y_pred += a(i);
        else if( i == 2 )
          y_pred += a(i) / x[bin];
        else
          y_pred += a(i) * std::pow( x[bin], static_cast<double>(i-1) );
      }
    }else
    {
      for( int i = 0; i < num_fit_coefficients; ++i )
        y_pred += a(i) * std::pow( x[bin], static_cast<double>(i) );
    }//if( include_inv_term ) / else
    
    y_pred = sqrt( y_pred );
    chi2 += std::pow( (y_pred - widths[bin]) / widths_uncert[bin], 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  if( include_inv_term && (coeffs.size() > 1) )
  {
    // Swap the zeroth and first coefficients so things are as expected.
    //  (i.e., constant term is zeroth, and energy dependant term is at index==1)
    std::swap( coeffs[0], coeffs[1] );
  }
  
  return chi2;
}//double fit_sqrt_poly_fwhm_lls(...)
  
  
double fit_intrinsic_eff_least_linear_squares( const std::vector<DetEffDataPoint> data,
                            const int order,
                           std::vector<float> &coeffs,
                           std::vector<float> &coeff_uncerts )
{
  const size_t nbin = data.size();
  
  //log(eff(x)) = A0 + A1*logx + A2*logx^2 + A3*logx^3
  vector<float> x, effs, effs_uncert;
  x.resize( data.size() );
  effs.resize( data.size() );
  effs_uncert.resize( data.size() );
  for( size_t i = 0; i < data.size(); ++i )
  {
    x[i] = data[i].energy;
    effs[i] = data[i].efficiency;
    effs_uncert[i] = data[i].efficiency_uncert;
  }
  
  
  //General Linear Least Squares fit
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation is quite inneficient
  //log(eff(x)) = A0 + A1*logx + A2*logx^2 + A3*logx^3 + ...
  using namespace boost::numeric;
  
  ublas::matrix<double> A( nbin, order );
  ublas::vector<double> b( nbin );
  
  for( size_t row = 0; row < nbin; ++row )
  {
    const double data_y = log(effs[row]);
    const double data_y_uncert = fabs( log(effs_uncert[row]) );
    
    b(row) = data_y / data_y_uncert;
    for( int col = 0; col < order; ++col )
      A(row,col) = std::pow( log(x[row]), double(col)) / data_y_uncert;
  }//for( int col = 0; col < order; ++col )
  
  const ublas::matrix<double> A_transpose = ublas::trans( A );
  const ublas::matrix<double> alpha = prod( A_transpose, A );
  ublas::matrix<double> C( alpha.size1(), alpha.size2() );
  const bool success = matrix_invert( alpha, C );
  if( !success )
    throw runtime_error( "fit_intrinsic_eff_least_linear_squares(...): trouble inverting matrix" );
  
  const ublas::vector<double> beta = prod( A_transpose, b );
  const ublas::vector<double> a = prod( C, beta );
  
  coeffs.resize( order );
  coeff_uncerts.resize( order );
  for( int coef = 0; coef < order; ++coef )
  {
    coeffs[coef] = static_cast<float>( a(coef) );
    coeff_uncerts[coef] = static_cast<float>( std::sqrt( C(coef,coef) ) );
  }//for( int coef = 0; coef < order; ++coef )
  
  double chi2 = 0;
  for( size_t bin = 0; bin < nbin; ++bin )
  {
    double y_pred = 0.0;
    for( int i = 0; i < order; ++i )
      y_pred += a(i) * std::pow( log(x[bin]), double(i) );
    y_pred = exp( y_pred );
    chi2 += std::pow( (y_pred - effs[bin]) / effs_uncert[bin], 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2;
}//double fit_intrinsic_eff_least_linear_squares(...)
  
  
double performEfficiencyFit( const std::vector<DetEffDataPoint> data,
                            const int fcnOrder,
                            std::vector<float> &result,
                            std::vector<float> &uncerts )
{
  if( data.empty() )
    throw runtime_error( "MakeDrfFit::performEfficiencyFit(...): no input peaks" );
  
  if( fcnOrder < 1 )
    throw runtime_error( "MakeDrfFit::performEfficiencyFit(...): requested fit order " + std::to_string(fcnOrder) + " (must be at least 1)" );
  
  if( fcnOrder > static_cast<int>(data.size()) )
    throw runtime_error( "MakeDrfFit::performEfficiencyFit(...): requested fit order " + std::to_string(fcnOrder)
                         + ", with only " + std::to_string(data.size()) + " data points." );
  
  result.clear();
  uncerts.clear();
  
  auto maxel = std::max_element( begin(data), end(data), [](const DetEffDataPoint &lhs,const DetEffDataPoint &rhs) -> bool {
    return lhs.energy > rhs.energy;
  } );
  
  DetectorEffFitness fitness( data, fcnOrder );
  ROOT::Minuit2::MnUserParameters inputPrams;
  
  
  const bool inMeV = (maxel->energy < 30);
  
  //Value limits taken from a wide variety of detectors, and then exaggerated
  //  to hopefully cover all reasonable ranges.
  const double mev_lower_bounds[8]    = {-4.75, -2.5,  -1.5,  -1.0,  -1.5,  -1.25, -0.75, -0.1 };
  const double mev_upper_bounds[8]    = {-0.75, -0.25,  1.25,  0.9,   0.25,  0.25,  0.05,  0.05 };
  const double mev_starting_values[8] = {-2.7,  -1.2,  -0.18, -0.14, -0.40, -0.15, -0.03, -0.002 };
  
  bool lls_worked = false;
  try
  {
    const double chi2 = fit_intrinsic_eff_least_linear_squares( data, fcnOrder, result, uncerts );
    
    cout << "GLLS Fit chi2=" << chi2 << ", coefs={";
    for( size_t i = 0; i < result.size(); ++i )
      cout << result[i] << "+-" << uncerts[i] << ", ";
    cout << "}" << endl;
    lls_worked = true;
    //return (fcnOrder == data.size()) ? chi2 : (chi2 / (data.size() - fcnOrder));
  
    //Least Linear Squares
    for( size_t i = 0; i < result.size(); ++i )
    {
      if( inMeV )
      {
        if( result[i] > mev_lower_bounds[i] && result[i] < mev_upper_bounds[i] )
          inputPrams.Add( std::to_string(i), result[i], (mev_upper_bounds[i] - mev_lower_bounds[i])/100.0, mev_lower_bounds[i], mev_upper_bounds[i] );
        else
          inputPrams.Add( std::to_string(i), result[i], (mev_upper_bounds[i] - mev_lower_bounds[i])/10.0 );
      }else
      {
        inputPrams.Add( std::to_string(i), result[i], 1.0 );
      }
    }
    
  }catch(std::exception &e)
  {
    cout << "GLLS Fit failed: " << e.what() << endl;
    
    if( inMeV )
    {
      for( int i = 0; i < 8 && i < fcnOrder; ++i )
      {
        assert( mev_upper_bounds[i] > mev_lower_bounds[i] );
        inputPrams.Add( std::to_string(i), mev_starting_values[i], (mev_upper_bounds[i] - mev_lower_bounds[i])/10.0, mev_lower_bounds[i], mev_upper_bounds[i] );
      }
      
      for( int order = 8; order <= fcnOrder; ++order )
        inputPrams.Add( std::to_string(order), 0, 0.001, -0.1, 0.1 );
    }else
    {
      inputPrams.Add( "A", -344, 100, -1000, 1000 );
      inputPrams.Add( "B", 270, 50 );
      inputPrams.Add( "C", -84, 50 );
      if( fcnOrder > 3 )
        inputPrams.Add( "D", 13.00, 10 );
      if( fcnOrder > 4 )
        inputPrams.Add( "E", -1.0, 0.1 );
      if( fcnOrder > 5 )
        inputPrams.Add( "F", 0.03, 0.01 );
      for( int order = 7; order <= fcnOrder; ++order )
        inputPrams.Add( std::string("")+char('A'+order-1), 0, 1.0 );
    }//if( inMeV ) / else
  }

  
  ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( fitness, inputParamState, strategy );
  
  double tolerance = 0.5;
  unsigned int maxFcnCall = 50000;
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, 2*tolerance );
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, static_cast<double>( data.size() ) );
  
  ROOT::Minuit2::MnUserParameters fitParams = fitter.Parameters();

  //Minos seems to give the same error we already have.
  //ROOT::Minuit2::MnMinos minos( fitness, minimum, 1 );
  //vector<pair<double,double>> assymerrors;
  //for( int i = 0; i < fcnOrder; ++i )
  //{
  //  ROOT::Minuit2::MinosError error = minos.Minos(i);
  //  assymerrors.push_back( make_pair(error.Lower(), error.Upper()) );
  //  cout << "Par " << i << " error: " << fitParams.Errors()[i] << ", Lower=" <<error.Lower() << ", Upper=" << error.Upper() << endl;
  //}
  
  
  double chi2 = fitness.DoEval( fitParams.Params() );
  const double lls_chi2 = lls_worked ? fitness.DoEval( vector<double>(begin(result),end(result)) ) : std::numeric_limits<double>::infinity();
  cerr << "Eff final chi2=" << chi2 << " vs " << lls_chi2 << " from LLS" << endl;
  
  if( !minimum.IsValid() )
  {
    stringstream msg;
    msg  << "Eff response function fit status is not valid"
    << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
    << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
    << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
    << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
    << "\n\tHasValidParameters: " << minimum.HasValidParameters()
    << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
    << endl;
    if( minimum.IsAboveMaxEdm() )
      msg << "\t\tEDM=" << minimum.Edm() << endl;
    cerr << endl << msg.str() << endl;
    
    if( !lls_worked )
      throw std::runtime_error( msg.str() );
  }//if( !minimum.IsValid() )
  
  if( lls_worked && (lls_chi2 < chi2) )
  {
    cerr << "Returning least linear squares answer for Intrinsic Eff equation" << endl;
    chi2 = lls_chi2;
  }else
  {
    result.clear();
    for( const double p : fitParams.Params() )
      result.push_back( static_cast<float>(p) );
  
    uncerts.clear();
    for( const double p : fitParams.Errors() )
      uncerts.push_back( p );
  }
  
  return (fcnOrder == data.size()) ? chi2 : (chi2 / (data.size() - fcnOrder));
}//performEfficiencyFit(...)
  
}//namespace MakeDrfFit
