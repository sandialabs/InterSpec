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

//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
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
        const float predicted = DetectorPeakResponse::peakResolutionSigma(
                                                                          peak->mean(), m_form, fx );
        if( predicted <= 0.0 || IsNan(predicted) || IsInf(predicted) )
          return 999999.9;
        
        chi2 += MakeDrfFit::peak_width_chi2( predicted, *peak );
      }//for( const EnergySigma &es : m_energies_and_sigmas )
      
      return chi2 / x.size();
    }//DoEval();
    
    virtual double operator()( const std::vector<double> &x ) const
    {
      return DoEval( x );
    }
    
    DetectorResolutionFitness&  operator=( const DetectorResolutionFitness &rhs )
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
}//namespace




namespace MakeDrfFit
{
  
double peak_width_chi2( double predicted_sigma, const PeakDef &peak )
{
  double measured_sigma = peak.sigma();
  double measured_sigma_uncert = peak.sigmaUncert();
  if( measured_sigma_uncert <= 0.01 )
    measured_sigma_uncert = 1.0;
  
  const double chi = (measured_sigma - predicted_sigma)/measured_sigma_uncert;
  return chi*chi;
}//peak_width_chi2(...)



void performResolutionFit( std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > peaks,
                           const size_t num_gamma_channels,
                           const DetectorPeakResponse::ResolutionFnctForm fnctnlForm,
                           std::vector<float> &answer,
                           std::vector<float> &uncerts )
{
  if( !peaks || peaks->empty() )
    throw runtime_error( "MakeDrfFit::performResolutionFit(...): no input peaks" );
  
  bool highres;
  if( num_gamma_channels )
  {
    highres = (num_gamma_channels > 3000);
  }else
  {
    const size_t index = peaks->size() / 2;
    const double ratio = (*peaks)[index]->fwhm() / (*peaks)[index]->mean();
    highres = (ratio < 0.03); //whatever, this is JIC anyway
  }//if( !!meas ) / else
  
  
  double a_initial, b_initial, c_initial, d_initial;
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
      
    case DetectorPeakResponse::kSqrtPolynomial:
    {
      //FWHM = A + B*pow( energy/1E3 + C*energy*energy/1E4, D );
      if( highres )
      {
        //Based on a single detector
        a_initial = 1.1;
        lowerA = -2.5;
        upperA = 5.0;
        
        b_initial = 0.7 / sqrt(0.661);
        lowerB = 0.0;
        upperB = 5.0/sqrt(0.661);
        
        c_initial = 0.04;
        lowerC = 0.0;
        upperC = 0.125;
      }else
      {
        //Based on zero detectors so far
        a_initial = 10;
        lowerA = -10;
        upperA = 20;
        
        b_initial = 60 / sqrt(0.661);
        lowerB = 0.0;
        upperB = 180/sqrt(0.661);
        
        c_initial = 0.04;
        lowerC = 0.0;
        upperC = 0.25;
      }//if( highres ) / else
      
      d_initial = 0.5;
      
      break;
    }//case kSqrtPolynomial:
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
      //We should never make it here if the detector FWHM functional form has
      //  not been explicitly set, so throw and error if we got here
      throw runtime_error( "MakeDrfFit::performResolutionFit(...):"
                          " invalid ResolutionFnctForm" );
      break;
  }//switch( fnctnlForm )
  
  answer.clear();
  uncerts.clear();
  
  //For high res we should only have to to this one fit
  //But for lowres detectors we should consider the cases of
  //  A==0, A<0, and A>0 when fitting for kGadrasResolutionFcn
  
  bool fitting_gad_a = false;
  DetectorResolutionFitness fitness( *peaks, fnctnlForm );
  
  ROOT::Minuit2::MnUserParameters inputPrams;
  switch( fnctnlForm )
  {
    case DetectorPeakResponse::kGadrasResolutionFcn:
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
      
    case DetectorPeakResponse::kSqrtPolynomial:
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
      break;
    }//case kSqrtPolynomial:
      
    case DetectorPeakResponse::kNumResolutionFnctForm:
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
  

  cerr << "FWHM final chi2=" << fitness.DoEval( fitParams.Params() ) << endl;
  
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
  
  for( size_t i = 0; i < fitParams.Params().size(); ++i )
    cout << "\tFit Par" << i << "=" << fitParams.Params()[i] << endl;
  
  for( const double p : fitParams.Params() )
    answer.push_back( static_cast<float>(p) );
  
  for( const double p : fitParams.Errors() )
    uncerts.push_back( p );
}//std::vector<float> performResolutionFit(...)

}//namespace MakeDrfFit
