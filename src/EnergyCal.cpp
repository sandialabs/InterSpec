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
#include <vector>
#include <numeric>
#include <cassert>
#include <sstream>
#include <iostream>
#include <stdexcept>

#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>


//Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinimize.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnUserParameterState.h"


#include "InterSpec/EnergyCal.h"
#include "SpecUtils/EnergyCalibration.h"

using namespace std;


namespace
{
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



/** As an alternative to #fit_energy_cal_poly (for use when, e.g., that fails)
 This class allows using Minuit to find energy calibration coefficients.
 This class may eventually go away after #fit_energy_cal_poly is tested enough.
 */
class PolyCalibCoefMinFcn
  : public ROOT::Minuit2::FCNBase
{
public:
  PolyCalibCoefMinFcn( const vector<EnergyCal::RecalPeakInfo> &peakInfo,
                       const size_t nbin,
                       SpecUtils::EnergyCalType eqnType,
                       const std::vector< std::pair<float,float> > &devpair )
  : ROOT::Minuit2::FCNBase(),
  m_nbin( nbin ),
  m_eqnType( eqnType ),
  m_peakInfo( peakInfo ),
  m_devpair( devpair )
  {
    switch( m_eqnType )
    {
      case SpecUtils::EnergyCalType::Polynomial:
      case SpecUtils::EnergyCalType::FullRangeFraction:
      case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
        break;
      case SpecUtils::EnergyCalType::LowerChannelEdge:
      case SpecUtils::EnergyCalType::InvalidEquationType:
        throw runtime_error( "PolyCalibCoefMinFcn can only wotrk with Full "
                            "Range Fraction and Polynomial binnings" );
        break;
    }
  }//PolyCalibCoefMinFcn( constructor )
  
  virtual ~PolyCalibCoefMinFcn(){}
  
  virtual double Up() const { return 1.0; }
  
  virtual double operator()( const std::vector<double> &coef ) const
  {
    double chi2 = 0.0;
    
    vector<float> float_coef;
    for( const double d : coef )
    {
      if( IsInf(d) || IsNan(d) )
      {
        fprintf( stderr, "Recalibrator::PolyCalibCoefMinFcn::operator() found invalid input parameter\n" );
        return 99999999.0;
      }
      
      float_coef.push_back( static_cast<float>(d) );
    }//for( const double d : coef )
    
    //  if( m_eqnType == SpecUtils::EnergyCalType::FullRangeFraction )
    //    float_coef = fullrangefraction_coef_to_polynomial( float_coef, m_nbin );
    if( m_eqnType == SpecUtils::EnergyCalType::Polynomial
       || m_eqnType == SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial )
      float_coef = SpecUtils::polynomial_coef_to_fullrangefraction( float_coef, m_nbin );
    
    
    for( const float d : float_coef )
    {
      if( IsInf(d) || IsNan(d) )
      {
        fprintf( stderr, "Recalibrator::PolyCalibCoefMinFcn::operator(): invalid conversion from Poly to FWF.\n" );
        return 99999999.0;
      }
    }//for( float d : float_coef )
    
    const float nearend = SpecUtils::fullrangefraction_energy( m_nbin-2, float_coef, m_nbin, m_devpair );
    const float end = SpecUtils::fullrangefraction_energy( m_nbin-1, float_coef, m_nbin, m_devpair );
    const float begin = SpecUtils::fullrangefraction_energy( 0, float_coef, m_nbin, m_devpair );
    const float nearbegin = SpecUtils::fullrangefraction_energy( 1, float_coef, m_nbin, m_devpair );
    
    if( (nearend >= end) || (begin >= nearbegin) )
    {
      fprintf( stderr, "Got nearend=%f, end=%f, begin=%f, nearbegin=%f\n", nearend, end, begin, nearbegin );
      return 99999999.0;
    }//if( almostLastEnergy > lastEnergy )
    
    for( const EnergyCal::RecalPeakInfo &info : m_peakInfo )
    {
      const double predictedMean = SpecUtils::fullrangefraction_energy( info.peakMeanBinNumber, float_coef, m_nbin, m_devpair );
      double uncert = ((info.peakMeanUncert<=0.0) ? 1.0 : info.peakMeanUncert );
      chi2 += pow(predictedMean - info.photopeakEnergy, 2.0 ) / (uncert*uncert);
    }//for( const &RecalPeakInfo info : peakInfo )
    
    if( IsInf(chi2) || IsNan(chi2) )
    {
      fprintf( stderr, "Recalibrator::PolyCalibCoefMinFcn::operator(): invalid result Chi2\n" );
      return 1000.0;
    }
    
    return chi2;
  }
  
protected:
  size_t m_nbin;
  SpecUtils::EnergyCalType m_eqnType;
  std::vector<EnergyCal::RecalPeakInfo> m_peakInfo;
  std::vector< std::pair<float,float> > m_devpair;
};//class PolyCalibCoefMinFcn

}//namespace


double EnergyCal::fit_energy_cal_poly( const std::vector<EnergyCal::RecalPeakInfo> &peakinfos,
                            const vector<bool> fitfor,
                            const std::vector<std::pair<float,float>> &dev_pairs,
                            vector<float> &coefs,
                            vector<float> &coefs_uncert )
{
  const size_t npeaks = peakinfos.size();
  const size_t fitorder = static_cast<size_t>( std::count(begin(fitfor),end(fitfor),true) );
  
  if( npeaks < 1 )
    throw runtime_error( "Must have at least one peak" );
  
  if( fitorder < 1 )
    throw runtime_error( "Must fit for at least one coefficient" );
  
  if( fitorder > npeaks )
    throw runtime_error( "Must have at least as many peaks as coefficients fitting for" );
  
  if( coefs.size() != fitfor.size() )
    throw runtime_error( "You must supply input coefficient when any of the coefficients are fixed" );
  
  //Energy = P0 + P1*x + P2*x^2 + P3*x^3, where x is bin number
  //  However, some of the coeffeicents may not be being fit for.
  vector<float> mean_bin( npeaks ), true_energies( npeaks ), energy_uncerts( npeaks );
  for( size_t i = 0; i < npeaks; ++i )
  {
    mean_bin[i] = peakinfos[i].peakMeanBinNumber;
    true_energies[i] = peakinfos[i].photopeakEnergy;
    energy_uncerts[i] = true_energies[i] * peakinfos[i].peakMeanUncert / std::max(peakinfos[i].peakMean,1.0);
  }
  
  //General Linear Least Squares fit
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation is quite inneficient
  //Energy_i = P0 + P1*pow(i,1) + P2*pow(i,2) + P3*pow(i,3)
  using namespace boost::numeric;
  
  ublas::matrix<double> A( npeaks, fitorder );
  ublas::vector<double> b( npeaks );
  
  for( size_t row = 0; row < npeaks; ++row )
  {
    double data_y = true_energies[row];
    const double data_y_uncert = fabs( energy_uncerts[row] );
    
    data_y -= SpecUtils::correction_due_to_dev_pairs( true_energies[row], dev_pairs );
    
    for( size_t col = 0, coef_index = 0; coef_index < fitfor.size(); ++coef_index )
    {
      if( fitfor[coef_index] )
      {
        assert( col < fitorder );
        A(row,col) = std::pow( mean_bin[row], double(coef_index)) / data_y_uncert;
        ++col;
      }else
      {
        data_y -= coefs[coef_index] * std::pow( mean_bin[row], double(coef_index));
      }
    }//
    
    b(row) = data_y / data_y_uncert;
  }//for( int col = 0; col < order; ++col )
  
  const ublas::matrix<double> A_transpose = ublas::trans( A );
  const ublas::matrix<double> alpha = prod( A_transpose, A );
  ublas::matrix<double> C( alpha.size1(), alpha.size2() );
  const bool success = matrix_invert( alpha, C );
  if( !success )
    throw runtime_error( "Trouble inverting least linear squares matrix" );
  
  const ublas::vector<double> beta = prod( A_transpose, b );
  const ublas::vector<double> a = prod( C, beta );
  
  coefs.resize( fitfor.size(), 0.0 );
  coefs_uncert.resize( fitfor.size(), 0.0 );
  
  for( size_t col = 0, coef_index = 0; coef_index < fitfor.size(); ++coef_index )
  {
    if( fitfor[coef_index] )
    {
      assert( col < fitorder );
      coefs[coef_index] = static_cast<float>( a(col) );
      coefs_uncert[coef_index] = static_cast<float>( std::sqrt( C(col,col) ) );
      ++col;
    }else
    {
      coefs_uncert[coef_index] = 0.0;
    }
  }//for( int coef = 0; coef < order; ++coef )
  
  double chi2 = 0;
  for( size_t bin = 0; bin < npeaks; ++bin )
  {
    double y_pred = 0.0;
    for( size_t i = 0; i < fitfor.size(); ++i )
      y_pred += coefs[i] * std::pow( mean_bin[bin], static_cast<double>(i) );
    y_pred += SpecUtils::deviation_pair_correction( y_pred, dev_pairs );
    chi2 += std::pow( (y_pred - true_energies[bin]) / energy_uncerts[bin], 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2;
}//double fit_energy_cal_poly(...)



double EnergyCal::fit_energy_cal_iterative( const std::vector<EnergyCal::RecalPeakInfo> &peakInfo,
                                           const size_t nbin,
                                           const SpecUtils::EnergyCalType eqnType,
                                           const std::vector<bool> fitfor,
                                           std::vector<float> &startingCoefs,
                                           const std::vector<std::pair<float,float>> &devpair,
                                           std::vector<float> &coefs,
                                           std::vector<float> &coefs_uncert,
                                           std::string &warning_msg )
{
  coefs.clear();
  coefs_uncert.clear();
  warning_msg.clear();
  
  if( peakInfo.size() < 1 )
    throw runtime_error( "fit_energy_cal_iterative: no peaks specified" );
  
  const size_t npars = startingCoefs.size();
  
  if( npars < 2 )
    throw runtime_error( "fit_energy_cal_iterative: must be at least coefficients" );
  
  if( fitfor.size() != npars )
    throw runtime_error( "fit_energy_cal_iterative: fitfor.size() != startingCoefs.size()" );
  
  switch( eqnType )
  {
    
    case SpecUtils::EnergyCalType::Polynomial:
    case SpecUtils::EnergyCalType::FullRangeFraction:
    case SpecUtils::EnergyCalType::UnspecifiedUsingDefaultPolynomial:
      break;
    case SpecUtils::EnergyCalType::LowerChannelEdge:
    case SpecUtils::EnergyCalType::InvalidEquationType:
      throw runtime_error( "fit_energy_cal_iterative: must be either Polynomial or"
                           " FullRangeFraction calibration type" );
      break;
  }//switch( eqnType )
  
  if( nbin < 7 )
    throw runtime_error( "fit_energy_cal_iterative: you must have at least 7 channels" );
  
  const size_t nfit = std::accumulate( begin(fitfor), end(fitfor), size_t(0) );
  //const size_t nfit = ([&]() -> size_t { size_t i = 0; for( auto f : fitfor ) i += f; return i; })();
  
  if( nfit < 1 )
    throw runtime_error( "fit_energy_cal_iterative: must fit for at least one coefficient" );
  
  PolyCalibCoefMinFcn chi2Fcn( peakInfo, nbin, eqnType, devpair );
  ROOT::Minuit2::MnUserParameters inputPrams;
    
  for( size_t i = 0; i < npars; ++i )
  {
    string name = std::to_string( i );
    double delta = 1.0/pow( static_cast<double>(nbin), static_cast<double>(i) );
    if( eqnType == SpecUtils::EnergyCalType::FullRangeFraction )
    {
      switch( i )
      {
        case 0: delta = 0.5; break;
        case 1: delta = 1.0; break;
        default:
          delta = 1.0 / std::pow( 10.0, static_cast<double>(i-1) );
          break;
      }//switch( i )
    }//if( FullRangeFraction )
        
    if( fitfor[i] )
      inputPrams.Add( name, startingCoefs[i], delta );
    else
      inputPrams.Add( name, startingCoefs[i] );
        
    if( i == 1 && fitfor[i] )
      inputPrams.SetLowerLimit( name, 0.0 );
  }//for( size_t i = 0; i < sm_numCoefs; ++i )
        
    
  ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
  ROOT::Minuit2::MnStrategy strategy( 1 ); //0 low, 1 medium, >=2 high
    
  const int nvarpars = inputPrams.VariableParameters();
  const unsigned int maxFcnCall = 200 + 100*nvarpars + 5*nvarpars*nvarpars;
  const double tolerance = 0.05; //0.05 * peakInfos.size();
    
  ROOT::Minuit2::CombinedMinimizer fitter;
  ROOT::Minuit2::FunctionMinimum minimum
    = fitter.Minimize( chi2Fcn, inputParamState,
                      strategy, maxFcnCall, tolerance );
    
  //Not sure why Minuit2 doesnt like converging on the minumum verry well, but
  //  rather than showing the user an error message, we'll give it anither try
  if( minimum.IsAboveMaxEdm() )
  {
    ROOT::Minuit2::MnMigrad fitter( chi2Fcn, inputParamState, strategy );
    minimum = fitter( maxFcnCall, tolerance );
  }//if( minimum.IsAboveMaxEdm() )
    
  if( !minimum.IsValid() )
  {
    if( !minimum.HasValidCovariance() || !minimum.HasValidParameters() )
    {
      string msg = "Fit for calibration parameters failed.";
      bool fithigher = false;
      for( size_t i = 2; i < npars; ++i )
        fithigher |= fitfor[i];
      
      if( fithigher )
        msg += " you might try not fitting for quadratic or higher terms.";
      
      throw runtime_error( msg );
    }//if( we are really messed up )
    
    cerr << minimum << endl;
    
    warning_msg = "Warning: calibration coefficient fit results "
    "may be invalid, please check, and if necessary "
    "revert this calibration instead.";
    
    if( minimum.IsAboveMaxEdm() )
    {
      warning_msg += " The estimated distance to optimal calibration parameters is "
      "too large.";
      
      cerr << "EDM=" << minimum.Edm() << ", tolerance=" << tolerance << ", npeaks=" << peakInfo.size() << endl;
    }
    if( minimum.HasReachedCallLimit() )
      warning_msg += " To many calls to chi2 routine were made.";
  }//if( !minimum.IsValid() )
  
  const ROOT::Minuit2::MnUserParameters &fitPrams = minimum.UserState().Parameters();
  const vector<double> parValues = fitPrams.Params();
  const vector<double> parErrors = fitPrams.Errors();
  assert( parValues.size() == npars );
  assert( parValues.size() == parErrors.size() );
  
  for( size_t i = 0; i < parValues.size(); ++i )
  {
    if( IsInf(parValues[i]) || IsNan(parValues[i]) )
      throw runtime_error( "Invalid calibration parameter from fit :(" );
  }
  
  string mmsg = "Minuit fit gave pars={";
  for( size_t i = 0; i < parValues.size(); ++i )
  {
    coefs.push_back( static_cast<float>(parValues[i]) );
    coefs_uncert.push_back( static_cast<float>(parErrors[i]) );
    
    mmsg += std::to_string(parValues[i]) + "+-" + std::to_string(parErrors[i]) + ", ";
  }//
  mmsg += "}\n\n";
  cout << mmsg;
  
  return minimum.Fval();
}//fit_energy_cal_iterative(...)
