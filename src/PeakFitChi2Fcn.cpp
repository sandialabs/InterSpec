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

// Block out some warnings occurring in boost files.
#pragma warning(disable:4800) // warning C4800: 'int' : forcing value to bool 'true' or 'false' (performance warning)

#include <memory>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/SpectrumChart.h"
#include "InterSpec/PeakFitChi2Fcn.h"

using namespace std;

static_assert( PeakDef::Mean == 0, "PeakDef::Mean != 0" );
static_assert( (PeakDef::Chi2DOF+1) == int(PeakDef::NumCoefficientTypes), "PeakDef::Chi2DOF not last enumerated value" );

static_assert( int(PeakFitChi2Fcn::Mean)             == int(PeakDef::Mean), "PeakFitChi2Fcn::Mean != PeakDef::Mean" );
static_assert( int(PeakFitChi2Fcn::Sigma)            == int(PeakDef::Sigma), "PeakFitChi2Fcn::Sigma != PeakDef::Sigma" );
static_assert( int(PeakFitChi2Fcn::GaussAmplitude)   == int(PeakDef::GaussAmplitude), "PeakFitChi2Fcn::GaussAmplitude != PeakDef::GaussAmplitude" );
static_assert( int(PeakFitChi2Fcn::SkewPar0)         == int(PeakDef::SkewPar0), "PeakFitChi2Fcn::SkewPar0 != PeakDef::SkewPar0" );
static_assert( int(PeakFitChi2Fcn::SkewPar1)         == int(PeakDef::SkewPar1), "PeakFitChi2Fcn::SkewPar1 != PeakDef::SkewPar1" );
static_assert( int(PeakFitChi2Fcn::SkewPar2)         == int(PeakDef::SkewPar2), "PeakFitChi2Fcn::SkewPar2 != PeakDef::SkewPar2" );
static_assert( int(PeakFitChi2Fcn::SkewPar3)         == int(PeakDef::SkewPar3), "PeakFitChi2Fcn::SkewPar3 != PeakDef::SkewPar3" );
static_assert( int(PeakFitChi2Fcn::Chi2DOF)          == int(PeakDef::Chi2DOF), "PeakFitChi2Fcn::Chi2DOF != PeakDef::Chi2DOF" );


std::atomic<size_t> MultiPeakFitChi2Fcn::sm_ncalls( 0 );
std::atomic<bool> MultiPeakFitChi2Fcn::sm_call_opt_integrate( false );

  
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
}


PeakFitChi2Fcn::PeakFitChi2Fcn( const int npeaks,
                                std::shared_ptr<const SpecUtils::Measurement> data,
                                std::shared_ptr<const SpecUtils::Measurement> continium )
  :  m_useReducedChi2( true ),
     m_useMultiPeakPunishment( false ),
     m_lower_channel(0),
     m_upper_channel(0),
     m_npeaks( npeaks ),
     m_data( data ),
     m_continium( continium )
{
}


PeakFitChi2Fcn::PeakFitChi2Fcn( const int npeaks,
                const size_t lower_channel, const size_t upper_channel,
                std::shared_ptr<const SpecUtils::Measurement> data,
                std::shared_ptr<const SpecUtils::Measurement> continium )
  :  m_useReducedChi2( true ),
     m_useMultiPeakPunishment( false ),
     m_lower_channel( lower_channel ),
     m_upper_channel( upper_channel ),
     m_npeaks( npeaks ),
     m_data( data ),
     m_continium( continium )
{
  const size_t nchannel = m_data ? m_data->num_gamma_channels() : size_t(0);
  if( nchannel < 3 )
    throw std::runtime_error( "PeakFitChi2Fcn: no spectrum" );
  
  if( m_lower_channel >= nchannel )
    m_lower_channel = nchannel;
  if( m_upper_channel >= nchannel )
    m_upper_channel = nchannel;
  
  if( m_lower_channel != m_upper_channel )
  {
    if( m_lower_channel >= nchannel )
      m_lower_channel = 0;
    if( m_upper_channel >= nchannel )
      m_upper_channel = nchannel - 1;
  }
}


double PeakFitChi2Fcn::Up() const
{
  return 1.0;
}

double PeakFitChi2Fcn::operator()( const double *x ) const
{
  return chi2( x );
}


double PeakFitChi2Fcn::operator()( const vector<double> &params ) const
{
  return operator()( &(params[0]) );
}


double PeakFitChi2Fcn::DoEval( const double *x ) const
{
  return operator()( x );
}


unsigned int PeakFitChi2Fcn::NDim() const
{
  return m_npeaks * NumFitPars;
}


PeakFitChi2Fcn *PeakFitChi2Fcn::Clone() const
{
  PeakFitChi2Fcn *drone = new PeakFitChi2Fcn( m_npeaks, m_data, m_continium );
  drone->m_useReducedChi2 = m_useReducedChi2;
  drone->m_useMultiPeakPunishment = m_useMultiPeakPunishment;
  drone->m_lower_channel = m_lower_channel;
  drone->m_upper_channel = m_upper_channel;

  return drone;
}//PeakFitChi2Fcn *Clone() const


PeakFitChi2Fcn &PeakFitChi2Fcn::operator=( const PeakFitChi2Fcn &rhs )
{
  if( this == &rhs )
    return *this;

  m_useReducedChi2 = rhs.m_useReducedChi2;
  m_useMultiPeakPunishment = rhs.m_useMultiPeakPunishment;
  m_lower_channel = rhs.m_lower_channel;
  m_upper_channel = rhs.m_upper_channel;
  m_npeaks = rhs.m_npeaks;
  m_data = rhs.m_data;
  m_continium = rhs.m_continium;

  return *this;
}



void PeakFitChi2Fcn::useReducedChi2( const bool use )
{
  m_useReducedChi2 = use;
}


void PeakFitChi2Fcn::setUseMultiPeakPunishment( bool punish )
{
  m_useMultiPeakPunishment = punish;
}


void PeakFitChi2Fcn::setSharedIndexToContinuumInfo( double &info, int index )
{
  if( index > 9998 )
    throw runtime_error( "Max number of peaks is 9998" );

  if( index < 0 )
    index = 9999;
  
  info = index + (10000.0 * std::floor(info / 10000.0));
}//setSharedIndexToContinuumInfo(...)


int PeakFitChi2Fcn::continuumInfoToSharedIndex( double info )
{
  info /= 10000.0;

  info = 10000.0 * (info - std::floor(info));
  
  //Note, the below "+ 0.5" is necessary, or else the cast to the int will
  //  occationally fail due to truncation due to impresise representation of
  //  the intermediate sub-integer values (e.g. when we divide by 10000.0).
  const int index = static_cast<int>( info + 0.5 );
  
  if( index == 9999 )
    return -1;
  return index;
}//continuumInfoToSharedIndex(...)

/*
bool PeakFitChi2Fcn::testOffsetConversions()
{
  for( PeakContinuum::OffsetType t = PeakContinuum::OffsetType(0);
       t <= PeakContinuum::External; t = PeakContinuum::OffsetType(t+1 ) )
  {
    for( int peakn = -1; peakn < 9999; peakn += 1 )
    {
      double info = 0.0;
      setSharedIndexToContinuumInfo( info, peakn );
      setOffsetTypeToContinuumInfo( info, t );
      int index = continuumInfoToSharedIndex( info );
      PeakContinuum::OffsetType type = continuumInfoToOffsetType( info );
      
      if( (type != t) || (index != peakn) )
        return false;
    }//for( int peakn = -1; peakn < 9999; peakn += 1 )
  }//for( PeakContinuum::OffsetType t...)

  //Now test setting values in the oppossit order
  for( PeakContinuum::OffsetType t = PeakContinuum::OffsetType(0);
      t <= PeakContinuum::External; t = PeakContinuum::OffsetType(t+1 ) )
  {
    for( int peakn = -1; peakn < 9999; peakn += 1 )
    {
      double info = 0.0;
      setOffsetTypeToContinuumInfo( info, t );
      setSharedIndexToContinuumInfo( info, peakn );
      int index = continuumInfoToSharedIndex( info );
      PeakContinuum::OffsetType type = continuumInfoToOffsetType( info );
      
      if( (type != t) || (index != peakn) )
        return false;
    }//for( int peakn = -1; peakn < 9999; peakn += 1 )
  }//for( PeakContinuum::OffsetType t...)
  
  return true;
}//void PeakFitChi2Fcn::testOffsetConversions()
*/

void PeakFitChi2Fcn::setOffsetTypeToContinuumInfo( double &info,
                                const PeakContinuum::OffsetType type )
{
  const double lower =  10000.0 * ((info/10000.0) - std::floor((info/10000.0)));
  const double upper = 100000.0 * floor(info/100000.0);
  const double digit = type * 10000.0;
  
  info = std::floor(lower+0.5) + digit + std::floor(upper+0.5);
}//setOffsetTypeToContinuumInfo(...)


PeakContinuum::OffsetType PeakFitChi2Fcn::continuumInfoToOffsetType( double info )
{
  info = std::floor( info / 10000.0 );
  info /= 10.0;
  info = 10.0 * (info - std::floor(info));
  
  return PeakContinuum::OffsetType( static_cast<int>(info) );
}//continuumInfoToOffsetType(...)


void PeakFitChi2Fcn::parametersToPeaks( std::vector<PeakDef> &peaks,
                                        const double *params,
                                        const double *errors ) const
{
  peaks.resize( m_npeaks );
  
  for( int peakn = 0; peakn < m_npeaks; ++peakn )
  {
    PeakDef &candidate_peak = peaks[peakn];
    const double *these_params = params + NumFitPars*peakn;
    const double *these_errors = errors ? errors + NumFitPars*peakn : 0;
  
    for( PeakDef::CoefficientType type = PeakDef::CoefficientType(0);
         type < PeakDef::NumCoefficientTypes;
         type = PeakDef::CoefficientType(type+1) )
    {
      candidate_peak.set_coefficient( these_params[type], type );
      if( these_errors )
        candidate_peak.set_uncertainty( these_errors[type], type );
    }
  
    candidate_peak.set_coefficient( 0.0, PeakDef::Chi2DOF );
    candidate_peak.set_uncertainty( 0.0, PeakDef::Chi2DOF );
    
    const double contInfo = these_params[ContinuumInfoField];
    const int sharedContIndex = continuumInfoToSharedIndex( contInfo );
    
    if( sharedContIndex >= 0 )
    {
      if( sharedContIndex >= peakn )
      {
        //shouldnt ever get here, this is more a test of programming logic
        cerr << "contInfo=" << contInfo
             << ", sharedContIndex=" << sharedContIndex << ", peakn=" << peakn
             << endl;
        for( FitPars p = FitPars(0); p < NumFitPars; p = FitPars(p+1) )
          cerr << "\t" << p << ": " << these_params[p] << endl;
        throw runtime_error( "PeakFitChi2Fcn::parametersToPeaks(...): invalid "
                             "shared continuum index" );
      }//if( sharedContIndex >= peakn )
      
      candidate_peak.setContinuum( peaks[sharedContIndex].continuum() );
    }else
    {
      PeakContinuum::OffsetType offset = continuumInfoToOffsetType( contInfo );
      const double lowx = these_params[RangeStartEnergy];
      const double highx = these_params[RangeEndEnergy];
      std::shared_ptr<PeakContinuum> continuum = candidate_peak.continuum();
      continuum->setRange( lowx, highx );
      continuum->setType( offset );
    
      if( continuum->isPolynomial() )
        continuum->setParameters( these_params[OffsetEnergyRelativeTo],
                         these_params + OffsetPolynomial0,
                         these_errors ? these_errors + OffsetPolynomial0 : 0 );
      else if( (offset==PeakContinuum::External) && m_continium )
        continuum->setExternalContinuum( m_continium );
    }//if( shares a continuum ) / else
    
    const double skew_type_d = round(these_params[SkewInfoField]);
    const auto skew_type_t = static_cast<PeakDef::SkewType>( static_cast<int>(skew_type_d) );
    bool valid_skew_type = false;
    switch( skew_type_t )
    {
      case PeakDef::NoSkew:
      case PeakDef::Bortel:
      case PeakDef::CrystalBall:
      case PeakDef::DoubleSidedCrystalBall:
      case PeakDef::GaussExp:
      case PeakDef::ExpGaussExp:
        valid_skew_type = true;
        break;
    }//switch( skew_type_t )
    
    assert( valid_skew_type );
    if( !valid_skew_type )
      throw std::logic_error( "Invalid skew type: " + std::to_string(these_params[SkewInfoField]) );
    
    candidate_peak.setSkewType( skew_type_t );
    
#if( PERFORM_DEVELOPER_CHECKS )
    if( m_data && (skew_type_t != PeakDef::NoSkew) )
    {
      const double mean = candidate_peak.mean();
      const double sigma = candidate_peak.sigma();
      
      const double lowx = these_params[RangeStartEnergy];
      const double highx = these_params[RangeEndEnergy];
      
      const size_t start_channel = m_data->find_gamma_channel(lowx); //m_data->find_gamma_channel(mean - 10*sigma);
      const size_t end_channel = m_data->find_gamma_channel(highx); //m_data->find_gamma_channel(mean + 10*sigma);
      assert( end_channel >= start_channel );
      const size_t nchannel = 1 + end_channel - start_channel;
      vector<double> counts( nchannel, 0.0 );
      auto energy_ptr = m_data->channel_energies();
      assert( energy_ptr );
      assert( energy_ptr->size() > end_channel );
      const float *energies = &(energy_ptr->at(start_channel));
      
      candidate_peak.gauss_integral(energies, &(counts[0]), nchannel );
      for( size_t i = 0; i < counts.size(); ++i )
      {
        if( IsInf(counts[i]) || IsNan(counts[i]) )
        {
          cerr << "Found invalid peak counts " << counts[i] << ", at " << energies[i] << " keV; parameter values:" << endl;
          
          for( PeakDef::CoefficientType ct = PeakDef::CoefficientType(0);
              ct < PeakDef::CoefficientType::NumCoefficientTypes;
              ct = PeakDef::CoefficientType(ct + 1) )
          {
            cerr << "\t" << std::setw(12) << PeakDef::to_string(ct) << ":" << candidate_peak.coefficient(ct) << endl;
          }
          cerr << endl;
        }
      }
    }//if( skew_type_t != PeakDef::NoSkew )
    
    //log_developer_error( __func__, "Invalid CSS color called back " );
#endif //PERFORM_DEVELOPER_CHECKS
  }//for( int peakn = 0; peakn < npeak; ++peakn )
}//parametersToPeak(...)


double PeakFitChi2Fcn::evalMultiPeakPunishment( const vector<const PeakDef *> &peaks ) const
{
  double chi2 = 0.0;
  const double punishment_chi2 = static_cast<double>( m_data->num_gamma_channels() );
  
  //punish if the peaks are too close.
  for( size_t i = 1; i < peaks.size(); ++i )
  {
    //we cant actually count on the peaks being sorted, so we just have to
    //  compare all the peaks (or rather all the preceeding peaks to catch all
    //  combinations)
    for( size_t j = 0; j < i; ++j )
    {
      const double sigma_i = peaks[i]->gausPeak() ? peaks[i]->sigma() : 0.25*peaks[i]->roiWidth();
      const double sigma_j = peaks[j]->gausPeak() ? peaks[j]->sigma() : 0.25*peaks[j]->roiWidth();
      
      const double sigma = 0.5*(sigma_j + sigma_i);
      const double dist = fabs(peaks[j]->mean() - peaks[i]->mean());
      double reldist = dist / sigma;
      if( reldist < 0.01 || IsInf(reldist) || IsNan(reldist) )
        reldist = 0.01;
      if( reldist < 1.0 )
        chi2 += (punishment_chi2 / reldist);
    }
  }//for( size_t i = 1; i < peaks.size(); ++i )
  
  //If the peak area is statistically insignificant on the interval
  //  -1.75sigma to 1.75, then punish!
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const double lower_energy = peaks[i]->gausPeak() ? peaks[i]->mean() - 1.75*peaks[i]->sigma() : peaks[i]->lowerX();
    const double upper_energy = peaks[i]->gausPeak() ? peaks[i]->mean() + 1.75*peaks[i]->sigma() : peaks[i]->upperX();
    
    //const size_t start_channel = m_data->find_gamma_channel( lower_energy );
    //const size_t end_channel = m_data->find_gamma_channel( upper_energy );
    //const double dataarea = m_data->gamma_channels_sum( lower_energy, upper_energy );
    const double dataarea = m_data->gamma_integral( lower_energy, upper_energy );
    
    if( peaks[i]->amplitude() < 2.0*sqrt(dataarea) )
    {
      if( peaks[i]->amplitude() <= 1.0 )
        chi2 += 100.0 * punishment_chi2;
      else
        chi2 += 0.5*(sqrt(dataarea)/peaks[i]->amplitude()) * punishment_chi2;
    }//if( peaks[i].amplitude() < 2.0*sqrt(dataarea) )
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  return chi2;
}//evalMultiPeakInsentive(...)


double PeakFitChi2Fcn::chi2( const double *params ) const
{
  assert( m_data );
  if( m_continium )
  {
    assert( m_data->num_gamma_channels() == m_continium->num_gamma_channels() );
  }

  const int nfitpar = m_npeaks * NumFitPars;
  for( int i = 0; i < nfitpar; ++i )
  {
    const double p = params[i];
    if( IsInf(p) || IsNan(p) )
    {
      static int ncalls = 0;
      if( ncalls++ < 10 )
        cerr << "\nPeakFitChi2Fcn::chi2(...): received invalid "
             << "input parameter " << i << endl;
      return DBL_MAX;
    }//if( IsInf(p) || IsNan(p) )
  }//for ( const double p : params )

  double chi2 = 0.0;
  int num_effective_bins = 0;
  
//  vector<double> gaussians( nbin, 0.0 );
  std::vector<PeakDef> peaks;
  parametersToPeaks( peaks, params );
  
  typedef map< std::shared_ptr<const PeakContinuum>, vector<const PeakDef *> > ContToPeakMap_t;
  ContToPeakMap_t contToPeakMap;
  
  for( int peakn = 0; peakn < m_npeaks; ++peakn )
  {
    const PeakDef &peak = peaks[peakn];
    std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
    if( continuum->type() == PeakContinuum::External )
      continuum.reset();
    contToPeakMap[continuum].push_back( &peaks[peakn] );
  }//for( size_t peak = 0; peak < m_npeaks; ++peak )
  
  for( const ContToPeakMap_t::value_type &vt : contToPeakMap )
  {
    const vector<const PeakDef *> &peaks = vt.second;
    const std::shared_ptr<const PeakContinuum> continuum = peaks[0]->continuum();
    
    std::set<size_t> binsToEval;
    vector<pair<size_t,size_t>> ranges; //the first channel number of range, and the number of channels in the range
    
    for( const PeakDef *peak : peaks )
    {
      size_t lower_channel, upper_channel;
      estimatePeakFitRange( *peak, m_data, lower_channel, upper_channel );
      
      if( m_lower_channel != m_upper_channel )
      {
        lower_channel = std::max( m_lower_channel, lower_channel );
        upper_channel = std::min( m_upper_channel, upper_channel );
      }//if( m_lowerbin != m_upperbin )
      
      for( size_t channel = lower_channel; channel <= upper_channel; ++channel )
        binsToEval.insert( channel );
      
      const size_t nchannel = 1 + (upper_channel - lower_channel);
      ranges.push_back( make_pair(lower_channel, nchannel) );
      
      if( continuum->energyRangeDefined() )
        break;
    }//for( const PeakDef *peak : peaks )
    
    if( (ranges.size() == 1) && m_data->channel_energies() )
    {
      const shared_ptr<const vector<float>> &channel_energies = m_data->channel_energies();
      assert( channel_energies );
      const float * const channel_energies_array = &((*channel_energies)[0]);
      
      const size_t nchannel = ranges[0].second;
      const size_t start_channel = ranges[0].first;
      const float * const energies = channel_energies_array + start_channel;
      
      vector<double> peak_counts( nchannel, 0.0 );
      
      for( const PeakDef *peak : peaks )
        peak->gauss_integral( energies, &(peak_counts[0]), nchannel );
      
      num_effective_bins += nchannel;
      for( size_t i = 0; i < nchannel; ++i )
      {
        const size_t channel = start_channel + i;
        const double ndata = m_data->gamma_channel_content(channel);
        const double nfitpeak = peak_counts[i];
        const double ncontinuum = continuum->offset_integral(energies[i], energies[i+1], m_data);
        
        if( ndata > 0.000001 )
          chi2 += pow( (ndata - ncontinuum - nfitpeak), 2.0 ) / ndata;
        else
          chi2 += fabs(nfitpeak + ncontinuum);  //This is a bit ad-hoc - is there a better solution? //XXX untested
      }
    }else
    {
      for( const size_t channel : binsToEval )
      {
        ++num_effective_bins;
        const double xbinlow = m_data->gamma_channel_lower(channel);
        const double xbinup = m_data->gamma_channel_upper(channel);
        const double ndata = m_data->gamma_channel_content(channel);
        const double ncontinuum = continuum->offset_integral(xbinlow, xbinup, m_data);
        
        double nfitpeak = 0.0;
        for( const PeakDef *peak : peaks )
          nfitpeak += peak->gauss_integral( xbinlow, xbinup );
        
        if( ndata > 0.000001 )
          chi2 += pow( (ndata - ncontinuum - nfitpeak), 2.0 ) / ndata;
        else
          chi2 += fabs(nfitpeak + ncontinuum);  //This is a bit ad-hoc - is there a better solution? //XXX untested
      }//for( int bin : binsToEval )
      
      if( m_useMultiPeakPunishment && (peaks.size() > 1) )
        chi2 += evalMultiPeakPunishment( peaks );
    }//for( const ContToPeakMap_t::value_type &vt : contToPeakMap )
  }//if( ranges.size() == 1 )
  
  if( m_useReducedChi2 )
  {
    // TODO: this is wrong number of fit parameters; doesnt account for parameters never fit, or parameters that are fixed
    const int nfitpar = m_npeaks * NumFitPars;
    const int ndof = ((num_effective_bins - nfitpar)>0)
                     ? (num_effective_bins - nfitpar)
                     : 1;
    chi2 /= ndof;
  }//if( m_useReducedChi2 )

  return chi2;
}//double chi2( const std::vector<double>& params ) const


void PeakFitChi2Fcn::addPeaksToFitter( ROOT::Minuit2::MnUserParameters &params,
                      const std::vector<PeakDef> &near_peaks,
                      std::shared_ptr<const SpecUtils::Measurement> data,
                      PeakFitChi2Fcn::AddPeaksToFitterMethod method )
{
  if( near_peaks.empty() )
    return;
  
  typedef map< std::shared_ptr<const PeakContinuum>, int > ContinuumMap_t;
  ContinuumMap_t continuumMap;
  
  double lowx, highx, fitrange;
  
  {
    double dummy;
    const PeakDef &lowgaus = near_peaks.front();
    const PeakDef &highgaus = near_peaks.back();
    findROIEnergyLimits( lowx, dummy, lowgaus, data );
    findROIEnergyLimits( dummy, highx, highgaus, data );
    fitrange = highx - lowx;
  }
  
  const int numStartPeaks = static_cast<int>((params.Params().size() / NumFitPars));
  
  for( size_t peaknum = 0; peaknum < near_peaks.size(); ++peaknum )
  {
    const PeakDef &peak = near_peaks[peaknum];
    if( peak.type() != PeakDef::GaussianDefined )
    {
      if( method == kFixPeakParameters )
        continue;
      
      throw runtime_error( "addPeaksToFitter cant add data defined peak to be fitted for" );
    }
    
    const int peakn = static_cast<int>( peaknum + numStartPeaks );
    
    if( params.Params().size() != peakn*NumFitPars )
      throw runtime_error( "Invalid logic adding parameters to fit" );
    
    const double sigma = peak.gausPeak() ? peak.sigma() : 0.25*peak.roiWidth();
    const double mean = peak.mean();
    const double amp = peak.amplitude();
    if( IsInf(mean) || IsNan(mean) || IsInf(sigma) || IsNan(sigma) || IsInf(amp)  || IsNan(amp) )
    {
      cerr << "\n\naddPeaksToFitter(): Peak with invalid mean, sigma, or amp"
           << endl;
      continue;
    }//if( problem with input peak )
    
    
    
    size_t lower_channel, upper_channel;
    estimatePeakFitRange( peak, data, lower_channel, upper_channel );
    const double lowerEdgeX = data->gamma_channel_lower(lower_channel);
    const double upperEdgeX = data->gamma_channel_upper(upper_channel);
    const double avrgbinwidth = (upperEdgeX - lowerEdgeX) / (1 + upper_channel - lower_channel);
    
    for( PeakDef::CoefficientType type = PeakDef::CoefficientType(0);
        type < PeakDef::NumCoefficientTypes;
        type = PeakDef::CoefficientType(type+1) )
    {
      const char *type_name = PeakDef::to_string( type );
      
      char name[128];
      snprintf( name, sizeof(name), "%s%i", type_name, peakn );
      
      double startval = 0.0, stepsize = 0.0, minval = 0.0, maxval = 0.0;
      
      switch( type )
      {
        case PeakDef::Mean:
          startval = peak.coefficient(type);
          stepsize = 0.1*peak.coefficient(type);
          
          if( method==kRefitPeakParameters )
          {
            minval = peak.coefficient(type) - 0.5*peak.coefficient(PeakDef::Sigma);
            maxval = peak.coefficient(type) + 0.5*peak.coefficient(PeakDef::Sigma);
          }else if( method == kFitUserIndicatedPeak )
          {
            minval = peak.mean() - 2.5*peak.sigma();
            maxval = peak.mean() + 2.5*peak.sigma();
            stepsize = 0.15 * peak.sigma();
            //should also so lower limits if there are nearby peaks
          }else
          {
            minval = peak.coefficient(type) - 2.0*peak.coefficient(PeakDef::Sigma);
            maxval = peak.coefficient(type) + 2.0*peak.coefficient(PeakDef::Sigma);
          }
          
          break;
          
        case PeakDef::Sigma:
          startval = peak.coefficient(type);
          stepsize = 0.1*peak.coefficient(type);
          
          if( method==kRefitPeakParameters )
          {
            minval = peak.coefficient(type) - 0.5*peak.coefficient(type);
            maxval = peak.coefficient(type) + 0.5*peak.coefficient(type);
          }else if( method == kFitUserIndicatedPeak )
          {
            const double sigma = peak.sigma();
            const double sigmaUncert = peak.sigmaUncert();
            
            if( sigmaUncert > 0.0 )
            {
              minval = 0.66*sigma;
              maxval = 1.5*sigma;
              stepsize = 0.05*sigma;
            }else
            {
              minval = 0.2*sigma;
              maxval = 2.5*sigma;
              stepsize = 0.25*sigma;
            }
          }else
          {
            minval = avrgbinwidth; //peak.coefficient(type) + avrgbinwidth;
            maxval = 0.45*fitrange; //peak.coefficient(type) - 0.45*fitrange
          }
          
          minval = std::max( minval, 0.15 );
          break;
          
        case PeakDef::GaussAmplitude:
          startval = peak.coefficient(type);
          stepsize = std::max( 0.1*peak.coefficient(type), sqrt(data->gamma_integral( mean-sigma, mean+sigma)) );
          minval = 0.0;
          maxval = data->gamma_count_sum();
          
          if( method == kFitUserIndicatedPeak )
            stepsize = 0.5*startval;
          break;
        
        case PeakDef::SkewPar0:
        case PeakDef::SkewPar1:
        case PeakDef::SkewPar2:
        case PeakDef::SkewPar3:
        {
          //We'll grab these later
       
          break;
        }//SkewPar...
          
        case PeakDef::Chi2DOF:
          startval = 0.0;
          stepsize = 0.0;
          break;
          
        case PeakDef::NumCoefficientTypes:
          break;
      }//switch( type )
      
      
      if( stepsize == 0.0 )
        stepsize = 0.1;
      
      if( method == kFixPeakParameters )
        startval = peak.coefficient(type);
      
      //The below commented-out code keeps from triggering an assert in Minuit2,
      //  however it doesnt take into account the parameters/situations we will never
      //  set minval/maxval for...
      //      if( (method != kFixPeakParameters) && (minval==maxval) )
      //      {
      //        string msg = "AddPeaksToFitter(): error determing extent of fit values "
      //                     "for fit type " + string(PeakDef::to_string(type));
      //        throw runtime_error( msg );
      //      }//if( Minuit2 will trigger an assert )
      
      switch( type )
      {
        case PeakDef::Mean:
        case PeakDef::Sigma:
        case PeakDef::GaussAmplitude:
          if( !peak.fitFor(type) )
          {
            params.Add( name,  startval );
          }else
          {
            switch( method )
            {
              case kFitForPeakParameters:
              case kFitUserIndicatedPeak:
                params.Add( name,  startval, stepsize, minval, maxval );
                break;
              case kFixPeakParameters:
                params.Add( name,  startval );
                break;
              case kRefitPeakParameters:
                params.Add( name,  startval, stepsize, minval, maxval );
                break;
            }//switch(method)
          }
          break;
          
          
        case PeakDef::SkewPar0:
        case PeakDef::SkewPar1:
        case PeakDef::SkewPar2:
        case PeakDef::SkewPar3:
        {
          bool do_fit = false;
          startval = 0.0, stepsize = 0.0, minval = 0.0, maxval = 0.0;
          const PeakDef::SkewType skew_type = peak.skewType();
          if( PeakDef::skew_parameter_range( skew_type, type, minval, maxval, startval, stepsize) )
          {
            do_fit = peak.fitFor(type);
            
            if( do_fit && (method == kRefitPeakParameters) )
            {
              if( (peak.coefficient(type) >= minval) && (peak.coefficient(type) <= maxval) )
              {
                startval = peak.coefficient(type);
                stepsize = std::max( 0.1*peak.coefficient(type), stepsize );
                minval = std::max( minval, fabs(0.5*peak.coefficient(type)) );
                maxval = std::min( maxval, 1.5*fabs(peak.coefficient(type)) );
              }//if( peaks current value is within the reasonable range )
              
              if( minval == maxval )
              {
                do_fit = false;
                startval = peak.coefficient(type);
              }
            }//if( method==kRefitPeakParameters )
          }//if( parameter replies )
          
          if( !do_fit || !peak.fitFor(type) )
          {
            params.Add( name,  peak.coefficient(type) );
          }else
          {
            switch( method )
            {
              case kFitForPeakParameters:
              case kFitUserIndicatedPeak:
              case kRefitPeakParameters:
                params.Add( name,  startval, stepsize, minval, maxval );
                break;
              case kFixPeakParameters:
                params.Add( name,  startval );
                break;
            }//switch(method)
          }//if( not fitting skew value ) / else
        }//case PeakDef::SkewPar...
          
        case PeakDef::Chi2DOF:
          params.Add( name, 0.0 );
          break;
          
        case PeakDef::NumCoefficientTypes:
          break;
      }//switch( type )
    }//for( loop over CoefficientType )
    
    int sharedContinuumIndex = -1;
    std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
    const ContinuumMap_t::iterator iter = continuumMap.find( continuum );
    if( iter != continuumMap.end() )
      sharedContinuumIndex = iter->second;
    else
      continuumMap[continuum] = static_cast<int>(peaknum + numStartPeaks);
    
    
    {// Begin set ContinuumInfoField
      double contInfo = 0.0;
      setSharedIndexToContinuumInfo( contInfo, sharedContinuumIndex );
      setOffsetTypeToContinuumInfo( contInfo, continuum->type() );
      
      const int resultIndex = continuumInfoToSharedIndex( contInfo );
      const PeakContinuum::OffsetType resultType
                                      = continuumInfoToOffsetType( contInfo );
      
      //Below check probably isn't necessary, but I'll leave in until tested
      //  or all platforms/architectures (implemented 20131230)
      if( resultIndex != sharedContinuumIndex
          || resultType != continuum->type() )
        throw runtime_error( "PeakFitChi2Fcn::addPeaksToFitter(...): logic "
                             "error encoding continuum info." );
      
      std::string name = "Peak" + std::to_string(peakn) + "ContInfo";
      params.Add( name,  contInfo );
    }// End set ContinuumInfoField
    
    
    {// Begin set SkewInfoField
      const double type = static_cast<int>( peak.skewType() );
      const string name = "Peak" + std::to_string(peakn) + "SkewType";
      params.Add( name,  type );
    }// End set SkewInfoField
    
    if( sharedContinuumIndex >= 0 )
    {
      string name = "StartEnergy" + std::to_string(peakn);
      params.Add( name, 0.0 );
      
      name = "EndEnergy" + std::to_string(peakn);
      params.Add( name, 0.0 );
      
      name = "EnergyRelTo" + std::to_string(peakn);
      params.Add( name,  0.0 );
      
      for( int i = 0; i < 5; ++i )
      {
        name = "Peak" + std::to_string(peakn)
               + "Par"+  std::to_string(i);
        params.Add( name,  0.0 );
      }
    }else
    {
      //add continuum coefficients to the fit
      //Prevois to 20131221, we would always use PeakDef::lowerX()/upperX(),
      //  where as now if this isnt defined, then the range will float in the
      //  fir, causing possible problems
      string name = "StartEnergy" + std::to_string(peakn);
      if( continuum->energyRangeDefined() )
        params.Add( name, peak.lowerX() );  //this will be the path taken for kFitUserIndicatedPeak, as of 20131230
      else
        params.Add( name, 0.0 );
    
      name = "EndEnergy" + std::to_string(peakn);
      if( continuum->energyRangeDefined() )
        params.Add( name, peak.upperX() );
      else
        params.Add( name, 0.0 );
      
      name = "EnergyRelTo" + std::to_string(peakn); //OffsetEnergyRelativeTo,
      params.Add( name,  continuum->referenceEnergy() );
      
      //We will always assume 4 polynomial coefficients
      const vector<double> polypars = continuum->parameters();
      const vector<bool> continuumFitForPar = continuum->fitForParameter();
      //    const PeakContinuum::OffsetType conttype = continuum->type();
      for( size_t i = 0; i < polypars.size(); ++i )
      {
        const double startval = polypars[i];
        double stepsize = 0.1*polypars[i];
        if( (method==kFitUserIndicatedPeak) || (stepsize <= 0.0) )
          stepsize = std::max( 0.25*startval, (100.0 / std::pow(10.0, 2.0*i)) );

        
        name = "Peak" + std::to_string(peakn) + "Par" + std::to_string(i);
      
        if( !continuumFitForPar[i] )
        {
          params.Add( name,  startval );
        }else
        {
          switch( method )
          {
            case kRefitPeakParameters:
            case kFitForPeakParameters:
            case kFitUserIndicatedPeak:
              params.Add( name,  startval, stepsize );
              //Keep the fitter form plunging off into the depths..., the
              //  continuum integral will never be allowed to be less 0.0,
              //  so to have a continuum, at least half of it should be
              //  above 0, or else why have one (note: subjective statment)
              if( i == 0 )
                params.SetLowerLimit( name, 0.0 );
              break;
            
            case kFixPeakParameters:
              params.Add( name,  startval );
              break;
          }//switch(method)
        }//if( peak.m_offsetType<PeakDef::Constant ) / else
      }//for( size_t i = 0; i < polypars.size(); ++i )
    
      for( size_t i = polypars.size(); i < 5; ++i )
      {
        name = "Peak" + std::to_string(peakn) + "Par" + std::to_string(i);
        params.Add( name,  0.0 );
      }
    }//if( shareContinuum ) / else
  }//for( add peaks to fitter )
}//void addPeaksToFitter( ROOT::Minuit2::Minuit2Minimizer &fitter, const std::vector<PeakDef> near_peaks )





/*
 Want to fit for a number of peaks in the given range.
 All peaks should share a common polynomial continuum.
 Widths of the peaks are tied together
 
 parameters:
 p0                        //continuum polynomial 0
 p1
 p2
 ...pX
 Skew0                     //Skew parameters, will be `PeakDef::num_skew_parameters(SkewType)`
 Skew1                     //  All peaks in the ROI share the same skew type/values.
 ...3
 width     1
 mean      1
 Amplitude 1
 width     2
 mean      2
 Amplitude 2
 ...
 width     X
 mean      X
 Amplitude X
 */
MultiPeakFitChi2Fcn::MultiPeakFitChi2Fcn( const int npeaks,
                                         std::shared_ptr<const SpecUtils::Measurement> data,
                                         PeakContinuum::OffsetType offsetType,
                                         PeakDef::SkewType skewType,
                                         const size_t lower_channel,
                                         const size_t upper_channel )
  : m_npeak( npeaks ),
    m_lower_channel( lower_channel ),
    m_upper_channel( upper_channel ),
    m_numOffset( 0 ),
    m_offsetType( offsetType ),
    m_skewType( skewType ),
    m_data( data ),
    m_reldiff_punish_start( 1.25 ),
    m_reldiff_punish_weight( 2.0 )
{
  if( !m_data )
    throw std::runtime_error( "MultiPeakFitChi2Fcn null data passed in" );
  
  const size_t nchannel = m_data->num_gamma_channels();
  if( nchannel < 3 )
    throw std::runtime_error( "MultiPeakFitChi2Fcn no gamma spectrum" );
  
  if( m_upper_channel >= nchannel )
    m_upper_channel = nchannel - 1;
  if( m_lower_channel >= nchannel )
    m_lower_channel = nchannel - 1;
  if( m_lower_channel > m_upper_channel )
    std::swap( m_lower_channel, m_upper_channel );
  
  m_rangeLow = data->gamma_channel_lower(m_lower_channel);
  m_highRange = data->gamma_channel_upper(m_upper_channel);
  
  m_nbin = 1 + m_upper_channel - m_lower_channel;
  m_energies.resize( m_nbin + 1 );
  m_dataCounts.resize( m_nbin );
  
  for( size_t channel = m_lower_channel; channel <= m_upper_channel; ++channel )
  {
    const size_t i = channel - m_lower_channel;
    m_energies[i] = m_data->gamma_channel_lower(channel);
    m_dataCounts[i] = m_data->gamma_channel_content(channel);
  }//for( int bin = m_lowerbin; bin <= m_upperbin; ++bin )
  
  assert( m_energies.size() == (m_nbin + 1) );
  m_energies[m_nbin] = m_data->gamma_channel_upper(m_upper_channel);
  
  switch( m_offsetType )
  {
    case PeakContinuum::NoOffset:
    case PeakContinuum::External:
      throw runtime_error( "MultiPeakFitChi2Fcn: invalid offset type" );
    break;
    
    case PeakContinuum::Constant:   case PeakContinuum::Linear:
    case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
      m_numOffset = static_cast<int>(m_offsetType);
    break;
    
    case PeakContinuum::FlatStep:
    case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
      m_numOffset = 2 + (m_offsetType - PeakContinuum::FlatStep);
    break;
  }//switch( m_offsetType )
  
  
  const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skewType );
  m_numOffset += static_cast<int>( num_skew_pars );
}//MultiPeakFitChi2Fcn constructor


MultiPeakFitChi2Fcn &MultiPeakFitChi2Fcn::operator=( const MultiPeakFitChi2Fcn &rhs )
{
  if( &rhs == this )
    return *this;
  
  m_npeak = rhs.m_npeak;
  m_lower_channel = rhs.m_lower_channel;
  m_upper_channel = rhs.m_upper_channel;
  m_numOffset = rhs.m_numOffset;
  m_nbin = rhs.m_nbin;
  m_rangeLow = rhs.m_rangeLow;
  m_highRange = rhs.m_highRange;
  m_energies = rhs.m_energies;
  m_dataCounts = rhs.m_dataCounts;
  m_offsetType = rhs.m_offsetType;
  m_skewType = rhs.m_skewType;
  m_data = rhs.m_data;
  m_reldiff_punish_start = rhs.m_reldiff_punish_start;
  m_reldiff_punish_weight = rhs.m_reldiff_punish_weight;
  
  return *this;
}//operator=


double MultiPeakFitChi2Fcn::Up() const
{
  return 1.0;
}


double MultiPeakFitChi2Fcn::operator()( const std::vector<double>& params ) const
{
  assert( params.size() == static_cast<size_t>(m_numOffset + 3*m_npeak) );
  
  if( sm_call_opt_integrate )
    return DoEval( &(params[0]), m_workingpeaks );
  return DoEval( &(params[0]) );
}



void MultiPeakFitChi2Fcn::parametersToPeaks( vector<PeakDef> &peaks,
                                             const double *x,
                                             const double *errors ) const
{
  peaks.resize( m_npeak );
  
  bool computeAreas = false;
  
  for( int peakn = 0; peakn < m_npeak; ++peakn )
  {
    PeakDef &peak = peaks[peakn];

    double sigma = x[m_numOffset + 3*peakn + 0];
    double centroid = x[m_numOffset + 3*peakn + 1];
    double amplitude = x[m_numOffset + 3*peakn + 2];
  
    if( peakn > 1 && sigma < -0.000001 )
      sigma = -sigma;
    else if( peakn > 1 && sigma > 0.000001 )
      throw runtime_error( "MultiPeakFitChi2Fcn::parametersToPeaks invalid width" );
    else if( m_npeak <= 1 )
      sigma = x[m_numOffset];
    else if( x[m_numOffset+3] > 0.000001 )
      sigma = x[m_numOffset] + ((centroid-m_rangeLow)/(m_highRange-m_rangeLow))*x[m_numOffset+3];
    else if( fabs(x[m_numOffset+3]) < 0.000001 )
      sigma = x[m_numOffset];
    else
      sigma = -x[m_numOffset + 3*peakn + 0];

    if( IsInf(sigma) || IsNan(sigma) )
      throw runtime_error( "MultiPeakFitChi2Fcn::parametersToPeaks: Invalid peak sigma." );
    
    if( IsInf(centroid) || IsNan(centroid) )
      throw runtime_error( "MultiPeakFitChi2Fcn::parametersToPeaks: Invalid peak centroid." );
    
    if( IsInf(amplitude) || IsNan(amplitude) )
      throw runtime_error( "MultiPeakFitChi2Fcn::parametersToPeaks: Invalid peak amplitude." );
    
    peak.setSigma( sigma );
    peak.setMean( centroid );
    peak.setAmplitude( amplitude );
   
    if( peakn )
    {
      peak.setContinuum( peaks[0].getContinuum() );
    }else
    {
      std::shared_ptr<PeakContinuum> continuum = peak.continuum();
      continuum->setType( m_offsetType );
      continuum->setRange( m_rangeLow, m_highRange );
      if( continuum->isPolynomial() )
        continuum->setParameters( m_rangeLow, x, errors );
    }//if( peakn ) / else
  
    if( fabs(amplitude+999.9) < 1.0 )
    {
      computeAreas = true;
      peak.setAmplitude( 1.0 );
    }else if( computeAreas )
    {
      throw runtime_error( "MultiPeakFitChi2Fcn: if you compute any areas from"
                           " data, you must compute all" );
    }//if( fabs(amplitude+999.9) < 0.5 ) / else
    
    if( errors )
    {
//      sigma = errors[m_numOffset + 3*peakn + 0];
      centroid = errors[m_numOffset + 3*peakn + 1];
      amplitude = errors[m_numOffset + 3*peakn + 2];
    
      if( peakn > 1 && sigma < -0.000001 )
        sigma = errors[m_numOffset + 3*peakn + 0];
      else if( m_npeak <= 1 )
        sigma = errors[m_numOffset];
      else if( x[m_numOffset+3] > 0.000001 )
        sigma = errors[m_numOffset] + ((centroid-m_rangeLow)/(m_highRange-m_rangeLow))*errors[m_numOffset+3];
      else if( fabs(x[m_numOffset+3]) < 0.000001 )
        sigma = errors[m_numOffset];
      else
        sigma = errors[m_numOffset + 3*peakn + 0];
      
      
    
      peak.setSigmaUncert( sigma );
      peak.setMeanUncert( centroid );
      peak.setAmplitudeUncert( amplitude );
    
      // \TODO: need to set continuum errors here
//      if( continuum->isPolynomial() )
    }//if( errors )
    
    peak.setSkewType( m_skewType );
    const size_t nskew = PeakDef::num_skew_parameters( m_skewType );
    const size_t ncont = PeakContinuum::num_parameters( m_offsetType );
    assert( static_cast<int>(nskew + ncont) == m_numOffset );
    
    for( size_t i = 0; i < nskew; ++i )
    {
      const double val = x[ncont + i];
      if( IsInf(val) || IsNan(val) )
        throw runtime_error( "MultiPeakFitChi2Fcn::parametersToPeaks: Invalid skew parameter value." );
      
      const auto par_type = static_cast<PeakDef::CoefficientType>( PeakDef::CoefficientType::SkewPar0 + i );
      peak.set_coefficient( val, par_type );
      if( errors )
        peak.set_uncertainty( errors[ncont + i], par_type );
    }//for( size_t i = 0; i < nskew; ++i )    
  }//for( int peakn = 0; peakn < m_npeak; ++peakn )
  
  if( computeAreas )
  {
    using namespace boost::numeric;
    ublas::matrix<double> f( m_npeak, m_npeak ), finv( m_npeak, m_npeak );
    ublas::vector<double> b( m_npeak ), a( m_npeak );
    
    vector<size_t> centroidbins( m_npeak, 0 );
    for( int i = 0; i < m_npeak; ++i )
    {
      size_t bin = lower_bound( m_energies.begin(), m_energies.end(),
                               static_cast<float>(peaks[i].mean()) ) - m_energies.begin();
      if( bin >= m_energies.size() )
        bin = m_energies.size() - 1;
      
      centroidbins[i] = bin;
      const double lowx = m_energies[bin];
      const double highx = m_energies[bin + 1];
      const double contarea = peaks[i].offset_integral( lowx, highx, m_data );
      
      b(i) = std::max( 0.0, m_dataCounts[bin] - contarea );
    }//for( int i = 0; i < m_npeak; ++i )
    
    for( int row = 0; row < m_npeak; ++row )
    {
      for( int col = 0; col < m_npeak; ++col )
      {
        const size_t bin = centroidbins[row];
        const double lowx = m_energies[bin];
        const double highx = m_energies[bin + 1];
        f(row,col) = peaks[col].gauss_integral( lowx, highx );
      }//for( int j = 0; j < peaks.size(); ++j )
    }//for( int i = 0; i < peaks.size(); ++i )
  
    bool success = false;
    try
    {
      success = matrix_invert( f, finv );
    }catch( std::exception &e )
    {
      cerr << "Failed matrix invert: " << e.what() << endl;
    }//try / catch
    
    if( success )
    {
      const ublas::vector<double> a = prod( finv, b );
      for( int i = 0; i < m_npeak; ++i )
      {
        if( IsInf(a(i)) || IsNan(a(i)) )
          throw runtime_error( "MultiPeakFitChi2Fcn::parametersToPeaks: Invalid computed peak area." );
        
        peaks[i].setAmplitude( a(i) );
//        cerr << "Setting candidate peak amplitude to " << a(i)
//             << " (would have set to " << (b(i) / f(i,i)) << ")" << endl;
      }
    }else
    {
      for( int i = 0; i < m_npeak; ++i )
      {
        double amp = (b(i) / f(i,i));
        amp = (IsInf(amp) || IsNan(amp)) ? 0.0 : amp;
        peaks[i].setAmplitude( amp );
      }
    }//if( success )
  }//if( computeAreas )
  
}//PeakDef parametersToPeaks( const double *x, int peakn )


double MultiPeakFitChi2Fcn::evalRelBinRange( const size_t beginRelChannel,
                                             const size_t endRelChannel,
                                             const std::vector<PeakDef> &peaks ) const
{
  if( peaks.empty() )
    throw runtime_error( "MultiPeakFitChi2Fcn::evalRelBinRange(): empty input" );
  
  assert( beginRelChannel < m_energies.size() );
  assert( endRelChannel <= m_energies.size() );
  
  if( (endRelChannel >= m_energies.size()) || (endRelChannel < beginRelChannel) )
    throw std::logic_error( "MultiPeakFitChi2Fcn::evalRelBinRange: invalid bin range" );
  
  double chi2 = 0.0;
  
  // We will get the peak contributions using the slightly more efficient methods
  const size_t nchan = endRelChannel - beginRelChannel;
  vector<double> peak_sum( endRelChannel - beginRelChannel, 0.0 );
  for( size_t i = 0; i < peaks.size(); ++i )
    peaks[i].gauss_integral( &(m_energies[beginRelChannel]), &(peak_sum[0]), nchan );
  
  for( size_t relchannel = beginRelChannel; relchannel < endRelChannel; ++relchannel )
  {
    double nfitpeak = peak_sum[relchannel - beginRelChannel];
    const double ndata = m_dataCounts[relchannel];
    const double xbinlow = m_energies[relchannel];
    const double xbinup  = m_energies[relchannel+1];
    const double ncontinuim = peaks[0].offset_integral( xbinlow, xbinup, m_data );
    
    const double datauncert = std::max( ndata, 1.0 );
    const double nabove = (ndata - ncontinuim - nfitpeak);
    chi2 += nabove*nabove / datauncert;
  }//for( int bin = xlowbin; bin <= xhighbin; ++bin )
  
  return chi2;
}//double MultiPeakFitChi2Fcn::evalRelBinRange


int MultiPeakFitChi2Fcn::nbin() const
{
  return static_cast<int>( m_nbin );
}

double MultiPeakFitChi2Fcn::dof() const
{
  const double num_continuum_pars = PeakContinuum::num_parameters( m_offsetType );
  
  // TODO: I'm not sure why its (num_continuum_pars - 1) below
  return 1.0*m_nbin - 2.0*m_npeak - (num_continuum_pars - 1);
}


void MultiPeakFitChi2Fcn::set_reldiff_punish_start( double reldiff_punish_start )
{
  m_reldiff_punish_start = reldiff_punish_start;
}//


void MultiPeakFitChi2Fcn::set_reldiff_punish_weight( double reldiff_punish_weight )
{
  m_reldiff_punish_weight = reldiff_punish_weight;
}//

double MultiPeakFitChi2Fcn::evalMultiPeakInsentive(
                                      const std::vector<PeakDef> &peaks ) const
{
  double chi2 = 0.0;
  const double punishment_chi2 = m_reldiff_punish_weight * m_nbin;
  
  //punish if the peaks are too close.
  for( size_t i = 1; i < peaks.size(); ++i )
  {
    for( size_t j = 0; j < i; ++j )
    {
      const double sigma = 0.5*(peaks[j].sigma() + peaks[i].sigma());
      const double dist = fabs(peaks[j].mean() - peaks[i].mean());
      double reldist = dist / sigma;
      if( reldist < 0.01 || IsInf(reldist) || IsNan(reldist) )
        reldist = 0.01;
      if( reldist < m_reldiff_punish_start )
        chi2 += (punishment_chi2 / reldist);
    }//for( size_t j = 0; j < i; ++j )
  }//for( size_t i = 1; i < peaks.size(); ++i )
  
  //If the peak area is statistically insignificant on the interval
  //  -1.75sigma to 1.75, then punish!
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    const double lower_energy = peaks[i].mean() - 1.75*peaks[i].sigma();
    const double upper_energy = peaks[i].mean() + 1.75*peaks[i].sigma();
    
//      const size_t peak_roi_begin = m_data->find_gamma_channel( lower_energy );
//      const size_t peak_roi_end = m_data->find_gamma_channel( upper_energy );
//      const double dataarea = m_data->gamma_channels_sum( peak_roi_begin, peak_roi_end );
    
    const size_t binstart = lower_bound( m_energies.begin(),
                                        m_energies.end(), static_cast<float>(lower_energy) )
                                         - m_energies.begin();
    const size_t binend = lower_bound( m_energies.begin(),
                                      m_energies.end(), static_cast<float>(upper_energy) )
                                       - m_energies.begin();
    double dataarea = 0.0;
    for( size_t bin = binstart; bin < binend; ++bin )
      dataarea += m_dataCounts[bin];
    
    if( peaks[i].amplitude() < 2.0*sqrt(dataarea) )
    {
      if( peaks[i].amplitude() <= 1 )
        chi2 += 100.0 * punishment_chi2;
      else
        chi2 += 0.5*(sqrt(dataarea)/peaks[i].amplitude()) * punishment_chi2;
    }//if( peaks[i].amplitude() < 2.0*sqrt(dataarea) )
    
//        const double sigma_ratio = fabs(peaks[i-1].sigma() - peaks[i].sigma())
//                                       / peaks[i].sigma();
//        const double sigma_dev = fabs(1.0-sigma_ratio);
//        if( sigma_dev > 0.05 )
//          chi2 += sigma_dev * punishment_chi2;
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  return chi2;
}//evalMultiPealInsentive(...)


double MultiPeakFitChi2Fcn::DoEval( const double *x ) const
{
  ++sm_ncalls;
  
  double chi2 = 0.0;
  vector<PeakDef> peaks( m_npeak );
  
  const size_t npars = static_cast<size_t>(m_numOffset + 3*m_npeak);
  for( size_t i = 0; i < npars; ++i )
    if( IsNan(x[i]) || IsInf(x[i]) )
      return DBL_MAX;
  
  parametersToPeaks( peaks, x );
  
  std::sort( peaks.begin(), peaks.end(), &PeakDef::lessThanByMean );
  
  chi2 += evalRelBinRange( 0, m_nbin, peaks );
  
  if( m_reldiff_punish_weight > 0.0 )
    chi2 += evalMultiPeakInsentive( peaks );
  
  if( IsInf(chi2) || IsNan(chi2) )
    return DBL_MAX;
  
  return chi2;
}//double DoEval( const double *x ) const


double MultiPeakFitChi2Fcn::DoEval( const double *x, std::vector<PeakDef> &peaks ) const
{
  ++sm_ncalls;
  
  double chi2 = 0.0;
  peaks.resize( m_npeak );
  
  const size_t npars = static_cast<size_t>(m_numOffset + 3*m_npeak);
  for( size_t i = 0; i < npars; ++i )
    if( IsNan(x[i]) || IsInf(x[i]) )
      return DBL_MAX;
  
  parametersToPeaks( peaks, x );
  
  std::sort( peaks.begin(), peaks.end(), &PeakDef::lessThanByMean );
  
  chi2 += evalRelBinRange( 0, m_nbin, peaks );
  
  if( m_reldiff_punish_weight > 0.0 )
    chi2 += evalMultiPeakInsentive( peaks );
  
  if( IsInf(chi2) || IsNan(chi2) )
    return DBL_MAX;
  
  return chi2;
}//double DoEval( const double *x ) const


LinearProblemSubSolveChi2Fcn::LinearProblemSubSolveChi2Fcn(
          const std::vector< std::shared_ptr<const PeakDef> > &originalPeaks,
          std::shared_ptr<const SpecUtils::Measurement> data,
          const PeakContinuum::OffsetType offsetType,
          const PeakDef::SkewType skewType,
          const float lowerROI, const float upperROI )
: ROOT::Minuit2::FCNBase(),
  m_nbin( 0 ),
  m_npeak( originalPeaks.size() ),
  m_lowerROI( lowerROI ),
  m_upperROI( upperROI ),
  m_offsetType( offsetType ),
  m_skewType( skewType ),
  m_originalPeaks( originalPeaks )
{
  init( data );
}


LinearProblemSubSolveChi2Fcn::LinearProblemSubSolveChi2Fcn( const size_t npeaks,
                               std::shared_ptr<const SpecUtils::Measurement> data,
                               const PeakContinuum::OffsetType offsetType,
                               const PeakDef::SkewType skewType,
                               const float lowerROI, const float upperROI )
: ROOT::Minuit2::FCNBase(),
  m_nbin( 0 ),
  m_npeak( npeaks ),
  m_lowerROI( lowerROI ),
  m_upperROI( upperROI ),
  m_offsetType( offsetType ),
  m_skewType( skewType ),
  m_originalPeaks{}
{
  init( data );
}

void LinearProblemSubSolveChi2Fcn::init( std::shared_ptr<const SpecUtils::Measurement> data )
{
  if( !data )
    throw runtime_error( "LinearProblemSubSolveChi2Fcn: invalid data" );
  
  if( m_npeak < 1 )
    throw runtime_error( "LinearProblemSubSolveChi2Fcn: invalid npeaks" );
  
  switch( m_offsetType )
  {
    case PeakContinuum::Constant:     case PeakContinuum::Linear:
    case PeakContinuum::Quadratic:    case PeakContinuum::Cubic:
    case PeakContinuum::FlatStep:     case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
    break;
      
    case PeakContinuum::NoOffset: case PeakContinuum::External:
      throw runtime_error( "LinearProblemSubSolveChi2Fcn: invalid offset" );
    break;
  }//switch( m_offsetType )
  
  const size_t lower_channel = data->find_gamma_channel( m_lowerROI );
  const size_t upper_channel = data->find_gamma_channel( m_upperROI );
  if( lower_channel > upper_channel )
    throw runtime_error( "LinearProblemSubSolveChi2Fcn: lower channel above upper channel" );
  
  m_nbin = 1 + upper_channel - lower_channel;
  
  if( m_nbin < 3 )
    throw runtime_error( "LinearProblemSubSolveChi2Fcn: invalid bin range" );
  
  m_y.reserve( m_nbin );
  m_x.reserve( m_nbin + 1 );
  for( size_t channel = lower_channel; channel <= upper_channel; ++channel )
  {
    m_y.push_back( data->gamma_channel_content(channel) );
    m_x.push_back( data->gamma_channel_lower(channel) );
  }//for( int bin = lowerbin; bin < upperbin; ++bin )
  
  //I dont think this next call should ever throw... if so its a programming error
  m_x.push_back( data->gamma_channel_lower(upper_channel+1) );
}//LinearProblemSubSolveChi2Fcn constructor



double LinearProblemSubSolveChi2Fcn::Up() const
{
  return 1.0;
}

size_t LinearProblemSubSolveChi2Fcn::nfitPars() const
{
  const size_t nskew = PeakDef::num_skew_parameters(m_skewType);
  if( m_npeak < 2 )
    return 2 + nskew;
  return m_npeak + 2 + nskew;
}//

double LinearProblemSubSolveChi2Fcn::operator()( const vector<double> &x ) const
{
  
  if( x.size() != nfitPars() )
    throw runtime_error( "LinearProblemSubSolveChi2Fcn::operator:"
                         " invalid number of parameters" );
  return DoEval( &x[0] );
}


double LinearProblemSubSolveChi2Fcn::DoEval( const double *x ) const
{
  double chi2 = DBL_MAX;
  try
  {
    std::vector<PeakDef> peaks;
    chi2 = parametersToPeaks( peaks, x );
    chi2 += punishment( peaks );
  }catch( std::exception &e )
  {
    cerr << "LinearProblemSubSolveChi2Fcn::DoEval caught: " << e.what() << endl;
    vector<double> means( m_npeak ), sigmas( m_npeak );
    for( size_t i = 0; i < m_npeak; ++i )
    {
      means[i] = x[i];
      const double frac = (x[i] - m_lowerROI) / (m_upperROI - m_lowerROI);
      sigmas[i] = x[m_npeak] + (m_npeak>1?(frac * x[m_npeak+1]) : 0.0);
    }//for( size_t i = 0; i < m_npeak; ++i )
    
    chi2 = 100.0*m_nbin;
    chi2 += 10.0*closenessPunishment( means, sigmas );
  }//try / catch
  
  return chi2;
}//DoEval(...)

size_t LinearProblemSubSolveChi2Fcn::nbin() const
{
  return m_nbin;
}

double LinearProblemSubSolveChi2Fcn::dof() const
{
  const double ncont_pars = PeakContinuum::num_parameters(m_offsetType);
  
  return 1.0*m_nbin - 2.0*m_npeak - ((ncont_pars > 0.0) ? (ncont_pars - 1.0) : 0.0);
}


double LinearProblemSubSolveChi2Fcn::parametersToPeaks( vector<PeakDef> &peaks,
                                                        const double *x,
                                                  const double *errors ) const
{
  peaks.clear();
  peaks.resize( m_npeak );
  for( size_t i = 1; i < m_npeak; ++i )
    peaks[i].setContinuum( peaks[0].continuum() );
  
  vector<size_t> indicesOfFittingPeaks;
  vector<PeakDef> fixedamppeaks;
  
  const size_t npar = nfitPars();
  const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skewType );
  const size_t start_skew_index = ((m_npeak < 2) ? size_t(2) : (2 + m_npeak));
  assert( start_skew_index <= npar );
  const double * const skew_params = x + start_skew_index;
  const double * const skew_param_errors = errors ? (errors + start_skew_index) : nullptr;
  
  vector<double> means, sigmas;
  for( size_t i = 0; i < m_npeak; ++i )
  {
    if( IsNan(x[i]) || IsInf(x[i]) )
      throw runtime_error( "NaN or Inf peak mean" );
    
    peaks[i].setMean( x[i] );
    if( errors )
      peaks[i].setMeanUncert( errors[i] );
    
    if( m_originalPeaks.empty()
       || m_originalPeaks[i]->fitFor(PeakDef::GaussAmplitude) )
    {
      means.push_back( x[i] );
    }
    
    peaks[i].setSkewType( m_skewType );
    for( size_t skew_par_index = 0; skew_par_index < num_skew_pars; ++skew_par_index )
    {
      const auto ct = PeakDef::CoefficientType(PeakDef::CoefficientType::SkewPar0 + skew_par_index);
      peaks[i].set_coefficient( skew_params[skew_par_index], ct );
      if( skew_param_errors )
        peaks[i].set_uncertainty( skew_param_errors[skew_par_index], ct );
    }//
  }//for( size_t i = 0; i < m_npeak; ++i )
  
  if( m_npeak < 2 )
  {
    double sigma = x[m_npeak];
    if( m_originalPeaks.size() && !m_originalPeaks[0]->fitFor(PeakDef::Sigma) )
      sigma = m_originalPeaks[0]->sigma();
    peaks[0].setSigma( sigma );
    if( errors )
      peaks[0].setSigmaUncert( errors[m_npeak] );
    
    if( m_originalPeaks.empty()
        || m_originalPeaks[0]->fitFor(PeakDef::GaussAmplitude) )
    {
      indicesOfFittingPeaks.push_back( 0 );
      sigmas.push_back( sigma );
    }else
    {
      peaks[0].setAmplitude( m_originalPeaks[0]->amplitude() );
      peaks[0].inheritUserSelectedOptions( *m_originalPeaks[0], true );
      fixedamppeaks.push_back( peaks[0] );
    }
  }else
  {
    for( size_t i = 0; i < m_npeak; ++i )
    {
      const double frac = (x[i] - m_lowerROI) / (m_upperROI - m_lowerROI);
      double sigma = x[m_npeak] + frac * x[m_npeak+1];
      
      if( m_originalPeaks.size() && !m_originalPeaks[i]->fitFor(PeakDef::Sigma) )
        sigma = m_originalPeaks[i]->sigma();
      
      peaks[i].setSigma( sigma );
      if( errors )
      {
        //XXX - I'm not entirely conviced uncert is correct (note: assuming
        //      100% correlation)
        const double uncert = errors[m_npeak] + frac * errors[m_npeak+1];
        peaks[i].setSigmaUncert( uncert );
      }//if( errors )
      
      if( m_originalPeaks.size()
         && !m_originalPeaks[i]->fitFor(PeakDef::GaussAmplitude) )
      {
        peaks[i].setAmplitude( m_originalPeaks[i]->amplitude() );
        peaks[i].setAmplitudeUncert( m_originalPeaks[i]->amplitudeUncert() );
        peaks[i].inheritUserSelectedOptions( *m_originalPeaks[i], true );
        fixedamppeaks.push_back( peaks[i] );
      }else
      {
        indicesOfFittingPeaks.push_back( i );
        sigmas.push_back( sigma );
      }
    }//for( size_t i = 0; i < m_npeak; ++i )
  }//if( one peak ) / else
  
  const int num_polynomial_terms = static_cast<int>( PeakContinuum::num_parameters(m_offsetType) );
  const bool step_continuum = PeakContinuum::is_step_continuum( m_offsetType );
  
  vector<double> amps, offsets, amps_uncerts, offsets_uncerts;
  
  const double chi2 = fit_amp_and_offset( &m_x[0], &m_y[0], m_nbin,
                                         num_polynomial_terms,
                                         step_continuum,
                                         m_lowerROI,
                                         means, sigmas,
                                         fixedamppeaks,
                                         m_skewType,
                                         skew_params,
                                         amps, offsets,
                                         amps_uncerts, offsets_uncerts );
  const double chi2Dof = chi2 / dof();
  
  peaks[0].continuum()->setType( m_offsetType );
  peaks[0].continuum()->setParameters( m_lowerROI, offsets, offsets_uncerts );
  peaks[0].continuum()->setRange( m_lowerROI, m_upperROI );
  
  for( size_t j = 0; j < indicesOfFittingPeaks.size(); ++j )
  {
    const size_t i = indicesOfFittingPeaks[j];
    peaks[i].setAmplitude( ((amps[j]>0.0) ? amps[j] : 0.0) );
    peaks[i].setAmplitudeUncert( amps_uncerts[j] );
  }//for( size_t i = 0; i < m_npeak; ++i )
  
  if( m_originalPeaks.size() )
  {
    for( size_t i = 0; i < m_npeak; ++i )
      peaks[i].inheritUserSelectedOptions( *m_originalPeaks[i], true );
  }//if( m_originalPeaks.size() )
  
  for( size_t j = 0; j < indicesOfFittingPeaks.size(); ++j )
  {
    const size_t i = indicesOfFittingPeaks[j];
    peaks[i].set_coefficient( chi2Dof, PeakDef::Chi2DOF );
  }
  
//  cerr << "\n\nFit:" << endl;
//  for( size_t i = 0; i < m_npeak; ++i )
//    cerr << peaks[i] << endl;
//  cerr << "\n\n" << endl;
  
  
  //More efficient to check if its sorted before sorting
  //  http://lexfridman.com/blogs/research/2012/07/04/sorting-a-sorted-list-in-c/
#if !ANDROID && __cplusplus > 199711L
  //is_sorted is c++11, apparently
  if( !std::is_sorted( peaks.begin(), peaks.end(), &PeakDef::lessThanByMean ) )
#endif
    std::sort( peaks.begin(), peaks.end(), &PeakDef::lessThanByMean );
  
  return chi2;
}//parametersToPeaks(...)

double LinearProblemSubSolveChi2Fcn::closenessPunishment(
                                          const vector<double> &means,
                                          const vector<double> &sigmas ) const
{
  double chi2 = 0.0;
  const double punishment_chi2 = 2.0*m_nbin;
  
  //punish if the peaks are too close.
  for( size_t i = 1; i < m_npeak; ++i )
  {
    for( size_t j = 0; j < i; ++j )
    {
      const double sigma = 0.5*(sigmas[j] + sigmas[i]);
      const double dist = fabs(means[j] - means[i]);
      double reldist = dist / sigma;
      if( reldist < 0.01 || IsInf(reldist) || IsNan(reldist) )
        reldist = 0.01;
      if( reldist < 1.25 )
        chi2 += (punishment_chi2 / reldist);
    }//for( size_t j = 0; j < i; ++j )
  }//for( size_t i = 1; i < peaks.size(); ++i )
  
  return chi2;
}//closenessPunishment(...)


void LinearProblemSubSolveChi2Fcn::addSkewParameters( ROOT::Minuit2::MnUserParameters &pars,
                              const PeakDef::SkewType skew_type )
{
  const size_t num_skew = PeakDef::num_skew_parameters( skew_type );
  for( size_t skew_num = 0; skew_num < num_skew; ++skew_num )
  {
    const PeakDef::CoefficientType ct = PeakDef::CoefficientType( PeakDef::SkewPar0 + skew_num );
    double lower_value, upper_value, starting_value, step_size;
    const bool use = PeakDef::skew_parameter_range( skew_type, ct, lower_value, upper_value,
                                                   starting_value, step_size );
    assert( use );
    if( !use )
      throw std::logic_error("inconsistent peak skew fit logic");
    
    const string name = "Skew" + std::to_string(skew_num);
    pars.Add( name, starting_value, step_size, lower_value, upper_value );
  }//for( size_t skew_num = 0; skew_num < num_skew; ++skew_num )
}//static void addSkewParameters(...);



PeakDef::SkewType LinearProblemSubSolveChi2Fcn::skewTypeFromPrevPeaks( const vector<shared_ptr<const PeakDef> > &inpeaks )
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
}//void skewFromPrevPeaks(...)


void LinearProblemSubSolveChi2Fcn::addSkewParameters( ROOT::Minuit2::MnUserParameters &pars,
                                                     const PeakDef::SkewType skew_type,
                                                     const vector<shared_ptr<const PeakDef> > &inpeaks )
{
  const size_t num_skew_pars = PeakDef::num_skew_parameters( skew_type );
  vector<bool> fit_parameter( num_skew_pars, true );
  vector<double> starting_value( num_skew_pars, 0.0 ), lower_values( num_skew_pars, 0.0 );
  vector<double> upper_values( num_skew_pars, 0.0 ), step_sizes( num_skew_pars, 0.0 );
  
  for( size_t i = 0; i < num_skew_pars; ++i )
  {
    const auto ct = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i);
    
    double lower, upper, start, dx;
    const bool use = PeakDef::skew_parameter_range( skew_type, ct, lower, upper, start, dx );
    assert( use );
    if( !use )
      throw logic_error( "Inconsistent skew par val" );
    
    starting_value[i] = start;
    lower_values[i] = lower;
    upper_values[i] = upper;
    step_sizes[i] = dx;
  }//for( size_t i = 0; i < num_skew_pars; ++i )
  
  
  for( const auto &p : inpeaks )
  {
    if( p->skewType() != skew_type )
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
      switch( skew_type )
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
      }//switch( skew_type )
      
      starting_value[i] = val;
    }//for( size_t i = 0; i < num_skew_pars; ++i )
  }//for( const auto &p : inpeaks )
  
  assert( num_skew_pars == fit_parameter.size() );
  assert( num_skew_pars == starting_value.size() );
  

  for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
  {
    const string name = "Skew" + std::to_string(skew_index);
    const double starting_val = starting_value[skew_index];
    
    if( fit_parameter[skew_index] )
    {
      const double lower = std::max( lower_values[skew_index], 0.5*starting_val );
      const double upper = std::min( upper_values[skew_index], 1.5*fabs(starting_val) );
      const double step = std::max( step_sizes[skew_index], 0.1*starting_val );
      
      pars.Add( name, starting_val, step, lower, upper );
    }else
    {
      pars.Add( name, starting_val );
    }
  }//for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
}//addSkewParameters(...)


double LinearProblemSubSolveChi2Fcn::punishment( const std::vector<PeakDef> &peaks ) const
{
  if( m_npeak != peaks.size() )
    throw runtime_error( "LinearProblemSubSolveChi2Fcn::punishment: invalid input" );
  
  vector<double> means( m_npeak ), sigmas( m_npeak ), amps( m_npeak );
  for( size_t i = 0; i < m_npeak; ++i )
  {
    means[i] = peaks[i].mean();
    sigmas[i] = peaks[i].sigma();
    amps[i] = peaks[i].amplitude();
  }//for( size_t i = 0; i < m_npeak; ++i )
  
  const double punishment_chi2 = 2.0*m_nbin;
  
  double chi2 = closenessPunishment( means, sigmas );
  
  //If the peak area is statistically insignificant on the interval
  //  -1.75sigma to 1.75, then punish!
  for( size_t i = 0; i < means.size(); ++i )
  {
    const double lower_energy = means[i] - 1.75*sigmas[i];
    const double upper_energy = means[i] + 1.75*sigmas[i];
        
    const size_t binstart = lower_bound( m_x.begin(), m_x.end(), lower_energy ) - m_x.begin();
    size_t binend = lower_bound( m_x.begin(), m_x.end(), upper_energy ) - m_x.begin();
	  binend = std::min( binend, m_y.size() );

    double dataarea = 0.0;
    for( size_t bin = binstart; bin < binend; ++bin )
      dataarea += m_y[bin];
    
    if( peaks[i].amplitude() < 2.0*sqrt(dataarea) )
    {
      if( amps[i] <= 1 )
        chi2 += 100.0 * punishment_chi2;
      else
        chi2 += 0.5*(sqrt(dataarea)/amps[i]) * punishment_chi2;
    }//if( peaks[i].amplitude() < 2.0*sqrt(dataarea) )
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  return chi2;
}//evalMultiPealInsentive(...)




