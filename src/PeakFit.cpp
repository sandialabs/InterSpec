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

#include <mutex>
#include <memory>
#include <limits>
#include <vector>
#include <utility>
#include <cstdlib>

#include <boost/math/constants/constants.hpp>

#define BOOST_UBLAS_TYPE_CHECK 0
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/math/distributions/poisson.hpp>


//Roots Minuit2 includes
#include "Minuit2/FCNBase.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnScan.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnSimplex.h"
#include "Minuit2/MinosError.h"
//#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/CombinedMinimizer.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/FumiliMinimizer.h"
#include "Minuit2/ScanMinimizer.h"
#include "Minuit2/SimplexMinimizer.h"
#include "Minuit2/MnMinimize.h"

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitChi2Fcn.h"
#include "SpecUtils/SpecUtilsAsync.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/DetectorPeakResponse.h"


#if( USE_REL_ACT_TOOL )
#include "InterSpec/PeakFitLM.h"
#endif

using namespace std;

using SpecUtils::Measurement;

// 20240911: The minimum uncertainty allowed for a gamma spectrum channel.
// Background subtracted spectra can end up with tiny bins, like 0.0007,
// which if we take its uncertainty to be its sqrt, a single bin like this will
// totally mess up the whole ROI.  So we'll impose a minimum uncertainty on
// each channel.
// However, if a spectrum is scaled, and not Poisson errored, this will totally
// mess things up (even though it wouldnt be right in the first place).
#define MIN_CHANNEL_UNCERT 1.0

template<class T>
bool matrix_invert( const boost::numeric::ublas::matrix<T>& input,
                   boost::numeric::ublas::matrix<T> &inverse )
{
  using namespace boost::numeric;
  ublas::matrix<T> A( input );
  ublas::permutation_matrix<std::size_t> pm( A.size1() );
  const size_t res = lu_factorize(A, pm);
  
  if( res != 0 )
  {
    //cout << "Singlular matrix passed in:" << endl;
    //for( size_t row = 0; row < input.size1(); ++row )
    //{
    //  for( size_t col = 0; col < input.size2(); ++col )
    //    cout << "\t" << input(row,col);
    //  cout << endl;
    //}
    return false;
  }//if( res != 0 )
  
  inverse.assign( ublas::identity_matrix<T>( A.size1() ) );
  lu_substitute(A, pm, inverse);
  return true;
}//matrix_invert


#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL )
namespace
{
  //inheriting from stringstream makes first string in line be printed as a
  //  pointer, not a string
  class DebugLog  //: public std::stringstream
  {
    //Make it so cout/cerr statments always end up non-interleaved when multiple
    // threads are calling cout/cerr
  public:
    explicit DebugLog( std::ostream &os ) : os(os) {}
    ~DebugLog() { os << ss.rdbuf() << std::flush; }
    template <typename T>
    DebugLog &operator<<(T const &t){ ss << t; return *this;}
    
  private:
    std::ostream &os;
    std::stringstream ss;
  };//class DebugLog
}//namespace
#endif

namespace ExperimentalAutomatedPeakSearch
{
  
bool largerByAmplitude( const std::shared_ptr<const PeakDef> &lhs, const std::shared_ptr<const PeakDef> &rhs )
{
  return lhs->amplitude() > rhs->amplitude();
}
  
  
void do_peak_automated_searchfit( const double x,
                                  const std::shared_ptr<const Measurement> &meas,
                                  const std::shared_ptr<const DetectorPeakResponse> &drf,
                                  const PeakShrdVec &inpeaks,
                                  const bool isHPGe,
                                  std::pair< PeakShrdVec, PeakShrdVec > &answer )
{
  try
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "Will try fitting peak clicked on at " << x << "\n";
#endif
    answer = searchForPeakFromUser( x, -1.0, meas, inpeaks, drf, isHPGe );
  }catch( std::exception &e )
  {
    cerr << "do_peak_searchfit(...): caught unexpected exception: '" << e.what()
         << "', returning blank results." << endl;
    answer.first.clear();
    answer.second.clear();
  }
}//void do_peak_searchfit(
  

  
std::vector<std::shared_ptr<const PeakDef> > filter_anomolous_width_peaks_highres(
                          const std::shared_ptr<const Measurement> meas,
                          std::vector<std::shared_ptr<const PeakDef> > input )
{
  if( !meas )
    return input;
  
  const size_t npeaks = input.size();
  
  if( npeaks < 3 )
    return input;
  
  if( npeaks <= 5 )
  {
    //We cant really fit for resolution here, but we can do some sanity checks
    //  To make sure peaks are thinner than higher energy peaks.
    double maxwidth = input.back()->sigma();
    double maxWidthMean = input.back()->mean();
    vector<std::shared_ptr<const PeakDef> > answer;
    
    for( int i = static_cast<int>(input.size()) - 1; i >= 0; --i )
    {
      const double sigma = input[i]->sigma();
      const double sigmaUncert = input[i]->sigmaUncert();
      
      if( sigma < 2.0*maxwidth )  //2.0 is made up
      {
        answer.push_back( input[i] );
        
        if( sigmaUncert > 0.0 && (sigmaUncert/sigma) < 0.2 && sigma > maxwidth )  //0.2 is completely made up
        {
          maxwidth = sigma;
          maxWidthMean = input[i]->mean();
        }
      }else
      {
        cerr << "\tRemoving peak at mean " << input[i]->mean()
             << ", sigma=" << sigma
             << " (for <=5 peak algo) for being wider than peak at "
             << maxWidthMean << " with sigma " << maxwidth << endl;
      }//if( sigma < 2.0*maxwidth ) / else
      
    }//for( int i = input.size()-1; i > 0; --i )
    
    
    return answer;
  }//if( npeaks <= 5 )
  
  
  DetectorPeakResponse det;
  
  std::shared_ptr<deque< std::shared_ptr<const PeakDef> > > inputpeaks( new std::deque< std::shared_ptr<const PeakDef> >() );
  for( size_t i = 0; i < npeaks; ++i )
  {
    if( fabs(511.0 - input[i]->mean()) > 10.0 )
      inputpeaks->push_back( input[i] );
  }
  
  try
  {
    det.fitResolution( inputpeaks, meas, DetectorPeakResponse::kSqrtPolynomial /*kGadrasResolutionFcn*/ );
  }catch( std::exception &e )
  {
    cerr << "filter_anomolous_width_peaks_highres: failed to fit resolution function: " << e.what() << endl;
    return input;
  }

  std::vector<std::shared_ptr<const PeakDef> > answer;
  
  for( size_t i = 0; i < npeaks; ++i )
  {
    if( fabs(511.0 - input[i]->mean()) < 10.0 )
    {
      answer.push_back( input[i] );
      continue;
    }
    
    const double mean = input[i]->mean();
    const double width = det.peakResolutionFWHM( mean );
    const double fracerror = fabs(width - input[i]->fwhm()) / std::min(input[i]->fwhm(), width);
    
    //For compton peak, could look ~2 FWHM on each side to see if decreasing...
    const bool tothin = ((width > input[i]->fwhm()) && fracerror > 0.8);
    const bool tothick = ((width < input[i]->fwhm()) && fracerror > 0.6);
    
    if( tothick || tothin )
    {
      //Actually we want to do some more sanity checks here....
      
      //If first or last, then fractional error should be 0.825
      
      //should actually try to fit 2 peaks in this region, and keep them if they
      //  have a better chi2
      
      int npeaksTouching = 0;
      for( size_t j = 0; j < npeaks; ++j )
        npeaksTouching += int(input[i]->continuum() == input[j]->continuum());
      
      if( npeaksTouching > 1 )
      {
        answer.push_back( input[i] );
        continue;
      }
      
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      cerr << "\tEliminating peak " << i << " of " << npeaks << " with mean="
      << mean << ", fractional error=" << fracerror
      << ", Actual FWHM=" << input[i]->fwhm() << ", fit resolution fwhm=" << width
      << endl;
#endif
      
      double twoPeaksChi2 = DBL_MAX;
      std::vector<std::shared_ptr<const PeakDef> > twopeaks;
      std::shared_ptr<const DetectorPeakResponse> detector;
      
      try
      {
        const double m = input[i]->mean();
        const double s = input[i]->sigma();

        PeakShrdVec onepeak( 1, input[i] );
        pair< PeakShrdVec, PeakShrdVec > twoPeaksPlus, twoPeaksMinus;
        twoPeaksPlus = searchForPeakFromUser( m + s, -1.0, meas, onepeak, nullptr, true );
        twoPeaksMinus = searchForPeakFromUser( m - s, -1.0, meas, onepeak, nullptr, true );
        
        if( twoPeaksPlus.first.size() == 2 && twoPeaksMinus.first.size() == 2 )
        {
          const double pChi2 = twoPeaksPlus.first[0]->chi2dof();
          const double mChi2 = twoPeaksMinus.first[0]->chi2dof();
          twopeaks = (pChi2 < mChi2) ? twoPeaksPlus.first : twoPeaksMinus.first;
        }else if( twoPeaksPlus.first.size() == 2 )
        {
          twopeaks = twoPeaksPlus.first;
        }else if( twoPeaksMinus.first.size() == 2 )
        {
          twopeaks = twoPeaksMinus.first;
        }
        
        if( twopeaks.size() != 2 )
          throw runtime_error( "Failed to fit two peaks in place of one" );
        
        twoPeaksChi2 = twopeaks[0]->chi2dof();
        
        if( twoPeaksChi2 <= (input[i]->chi2dof() + 1.2)
            && (twopeaks[0]->sigma() < input[i]->sigma()) )
        {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          cerr << "\tAdding two peaks in its place with chi2=" << twoPeaksChi2
               << " (single peak chi2=" << input[i]->chi2dof() << ")"
               << ", peak0(mean,sigma,amp)={" << twopeaks[0]->mean()
               << "," << twopeaks[0]->sigma() << "," << twopeaks[0]->amplitude()
               << "} (expected width " << det.peakResolutionSigma( twopeaks[0]->mean() )
               << "), and peak1(mean,sigma,amp)={" << twopeaks[1]->mean()
               << "," << twopeaks[1]->sigma() << "," << twopeaks[1]->amplitude()
               << "} (expected width " << det.peakResolutionSigma( twopeaks[1]->mean() )
               << ")"
               << endl;
#endif
          answer.push_back( twopeaks[0] );
          answer.push_back( twopeaks[1] );
        }//if( twoPeaksChi2 <= input[i]->chi2dof() )
        
        throw runtime_error( "Chi2 didnt improve" );
      }catch( std::exception &e )
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        cerr << "Failed to put two peaks where there was a wide one:"
             << e.what() << endl;
#endif
        twoPeaksChi2 = DBL_MAX;
      }//try / catch
    }else
    {
      // cerr << "\tPeak at " << mean << " has fractional error " << fracerror
      // << ", FWHM=" << input[i]->fwhm() << " fit peak width fwhm=" << width
      // << endl;
      answer.push_back( input[i] );
    }
  }//for( size_t i = 0; i < npeaks; ++i )

  
  return answer;
}//filter_anomolous_width_peaks_highres(...)

  
  

std::vector<std::shared_ptr<const PeakDef> > search_for_peaks_multithread(
                                       const std::shared_ptr<const Measurement> meas,
                                       const std::shared_ptr<const DetectorPeakResponse> &drf,
                                       std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > origpeaks,
                                       const bool isHPGe )
{
  typedef std::shared_ptr<PeakDef> PeakPtr;
  typedef std::shared_ptr<const PeakDef> PeakConstPtr;
  
  size_t lower_channel = 0, upper_channel = 0;
  //    ExperimentalPeakSearch::find_spectroscopic_extent( meas, lower_channel, upper_channel );
  //    cout << "Start at " << meas->gamma_channel_center( lower_channel ) << " and going through "
  //    << meas->gamma_channel_center( upper_channel ) << endl;
  
  const vector<std::shared_ptr<PeakDef> > initialcandidates
    = secondDerivativePeakCanidatesWithROI( meas, isHPGe, lower_channel, upper_channel );
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  {
    DebugLog log(cout);
    log << "multithread Found means are { ";
    for( size_t i = 0; i < initialcandidates.size(); ++i )
      log << (i?", ":"") << initialcandidates[i]->mean();
    log << " }\n";
  }
#endif
  
  vector<std::shared_ptr<const PeakDef> > fitpeakvec;
  
  if( !!origpeaks )
  {
    for( const PeakConstPtr &p : *origpeaks )
      fitpeakvec.push_back( p );
  }//if( !!origpeaks )
  
  vector<std::shared_ptr<const PeakDef> > candidates;
  if( !fitpeakvec.empty() )
  {
    for( const PeakPtr &p : initialcandidates )
    {
      bool isnear = false;
      for( const PeakConstPtr &orig : fitpeakvec )
      {
        if( orig->gausPeak() )
        {
          const double meandiff = orig->mean() - p->mean();
          const double minsigma = std::min( orig->sigma(), p->sigma() );
          isnear |= (fabs(meandiff/minsigma) < 0.75);
        }
      }//for( const PeakConstPtr &orig : fitpeakvec )
      
      if( !isnear )
        candidates.push_back( p );
    }//for( const PeakDef &p : candidates )
  }else
  {
    for( const PeakPtr &p : initialcandidates )
      candidates.push_back( p );
  }//if( !fitpeakvec.empty() ) / else

  
  while( !candidates.empty() )
  {
    std::sort( candidates.begin(), candidates.end(), &largerByAmplitude );
    
    //So this is kinda messy.  If we select a candidate peak for fitting, we
    //  have to not only make sure that neighboring candidate peaks wont also be
    //  fit at in this same iteration of the while loop (in another thread), but
    //  we have to make sure the peaks that the potentially connected candidate
    //  peaks are connected to, wont be fit - keeping in mind peaks in
    //  'fitpeakvec' can cause two previously un-connected candidate peaks, to
    //  now become connected.
    vector<std::shared_ptr<const PeakDef> > candidatesBeingFitFor;
    vector<std::shared_ptr<const PeakDef> > candidatesNotFitForDueToCausality;
    
    
    for( size_t i = 0; i < candidates.size(); ++i )
    {
      std::shared_ptr<const PeakDef> peak = candidates[i];
      
      if( std::find( candidatesNotFitForDueToCausality.begin(),
                     candidatesNotFitForDueToCausality.end(), peak )
          != candidatesNotFitForDueToCausality.end() )
        continue;
      
      vector<std::shared_ptr<const PeakDef> > inpeaks;
//      inpeaks.insert( inpeaks.end(), candidatesBeingFitFor.begin(), candidatesBeingFitFor.end() );
      inpeaks.insert( inpeaks.end(), candidatesNotFitForDueToCausality.begin(), candidatesNotFitForDueToCausality.end() );
      inpeaks.insert( inpeaks.end(), fitpeakvec.begin(), fitpeakvec.end() );
      inpeaks.insert( inpeaks.end(), candidates.begin()+i+1, candidates.end() );
      inpeaks.push_back( peak );
      
      //const double nsigma = isHPGe ? 10.0 : 5.0;
      const double nsigma = 7.5;
      const vector< vector<std::shared_ptr<const PeakDef> > > disconnectedpeaks
                              = causilyDisconnectedPeaks(  nsigma, true, inpeaks );
      
      for( size_t j = 0; j < disconnectedpeaks.size(); ++j )
      {
        const vector<std::shared_ptr<const PeakDef> > &peaks = disconnectedpeaks[j];
        if( std::find( peaks.begin(), peaks.end(), peak ) == peaks.end() )
          continue;
        
        for( const PeakConstPtr &p : peaks )
        {
          if( p == peak )
            continue;
          
          if( std::find( fitpeakvec.begin(), fitpeakvec.end(), p ) != fitpeakvec.end() )
            continue;
          
          if( std::find( candidatesNotFitForDueToCausality.begin(),
                         candidatesNotFitForDueToCausality.end(), p )
              != candidatesNotFitForDueToCausality.end() )
            continue;
          
          candidatesNotFitForDueToCausality.push_back( p );
        }
      }//for( size_t j = 0; j < disconnectedpeaks.size(); ++j )
      
      candidatesBeingFitFor.push_back( peak );
    }//for( size_t i = 0; i < candidates.size(); ++i )
    
    SpecUtilsAsync::ThreadPool pool;
    
    vector< pair< PeakShrdVec, PeakShrdVec > > results( candidatesBeingFitFor.size() );
    
    for( size_t i = 0; i < candidatesBeingFitFor.size(); ++i )
    {
      const double mean = candidatesBeingFitFor[i]->mean();
      pool.post( boost::bind( &do_peak_automated_searchfit, mean,
                             boost::cref(meas), boost::cref(drf),
                             boost::cref(fitpeakvec), 
                             isHPGe,
                             boost::ref(results[i]) ));
    }//for( size_t i = 0; i < peaksToTryIndices.size(); ++i )
    
    pool.join();
    
    bool was_collision = false;
    
    for( size_t i = 0; i < results.size(); ++i )
    {
      const PeakShrdVec &toadd = results[i].first;
      const PeakShrdVec &toremove = results[i].second;

      for( const PeakConstPtr &p : toremove )
      {
        vector<std::shared_ptr<const PeakDef> >::iterator pos
                            = std::find(fitpeakvec.begin(),fitpeakvec.end(),p);
        if( pos != fitpeakvec.end() )
          fitpeakvec.erase( pos );
        else
        {
          cerr << "There was a collision in the fit" << endl;
          was_collision = true;
        }
      }
      
      for( const PeakConstPtr &p : toadd )
        fitpeakvec.push_back( p );
    }//for( size_t i = 0; i < peaksToTryIndices.size(); ++i )
  
    std::sort( fitpeakvec.begin(), fitpeakvec.end(),
            &PeakDef::lessThanByMeanShrdPtr );
  
    vector<PeakConstPtr> nextcandidates;
    for( const PeakConstPtr &peak : candidates )
    {
      if( !std::count(candidatesBeingFitFor.begin(), candidatesBeingFitFor.end(), peak) )
        nextcandidates.push_back( peak );
    }//for( size_t i = 0; i < candidates.size(); ++i )
  
    candidates.swap( nextcandidates );
  }//while( !candidates.empty() )
  
  const auto detResolution = PeakFitUtils::coarse_resolution_from_peaks( fitpeakvec );
  if( detResolution == PeakFitUtils::CoarseResolutionType::High )
    fitpeakvec = filter_anomolous_width_peaks_highres( meas, fitpeakvec );
  
  return fitpeakvec;
}//search_for_peaks_multithread(...)
  
  

vector<std::shared_ptr<const PeakDef> > search_for_peaks_singlethread(
                        const std::shared_ptr<const Measurement> meas,
                        const std::shared_ptr<const DetectorPeakResponse> &drf,
                        std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > origpeaks,
                        const bool isHPGe )
{
  typedef std::shared_ptr<PeakDef> PeakPtr;
  typedef std::shared_ptr<const PeakDef> PeakConstPtr;
  
  size_t lower_channel = 0, upper_channel = 0;
  vector<PeakPtr> candidates
   = secondDerivativePeakCanidatesWithROI( meas, isHPGe, lower_channel, upper_channel );
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  {
    DebugLog log(cout);
    log << "Found means are { ";
    for( size_t i = 0; i < candidates.size(); ++i )
      log << (i?", ":"") << candidates[i]->mean();
    log << " }\n";
  }
#endif
  
  vector<std::shared_ptr<const PeakDef> > fitpeakvec;
  
  if( !!origpeaks )
  {
    for( const PeakConstPtr &p : *origpeaks )
      fitpeakvec.push_back( p );
  }//if( !!origpeaks )
  
  
  while( !candidates.empty() )
  {
    size_t largest_index = 0;
    for( size_t i = 1; i < candidates.size(); ++i )
    {
      if( candidates[i]->amplitude() > candidates[largest_index]->amplitude() )
        largest_index = i;
    }//for( size_t i = 1; i < candidates.size(); ++i )
    
    const PeakDef p = *candidates[largest_index];
    candidates.erase( candidates.begin() + largest_index );
    
    if( !!origpeaks )
    {
      bool originalnear = false;
      for( const PeakConstPtr &orig : *origpeaks )
        originalnear |= (orig->gausPeak() && fabs((orig->mean()-p.mean())/std::min(orig->sigma(),p.sigma())) < 0.75);
      
      if( originalnear )
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cout) << "Not trying to click on candidate peak at " << p.mean()
        << " keV since there was originally a peak there already\n";
#endif
        continue;
      }
    }//if( !!origpeaks )
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "Will try clicking on " << p.mean() << "\n";
#endif
    
    pair< PeakShrdVec, PeakShrdVec > results;
    do_peak_automated_searchfit( p.mean(), meas, drf, fitpeakvec, isHPGe, results );
    
    const PeakShrdVec &toadd = results.first;
    const PeakShrdVec &toremove = results.second;
    
    
    PeakShrdVec::iterator pos;
    for( size_t i = 0; i < toremove.size(); ++i )
    {
      pos = std::find( fitpeakvec.begin(), fitpeakvec.end(), toremove[i] );
      if( pos != fitpeakvec.end() )
        fitpeakvec.erase( pos );
    }
    
    for( size_t i = 0; i < toadd.size(); ++i )
      fitpeakvec.push_back( toadd[i] );
    
    std::sort( fitpeakvec.begin(), fitpeakvec.end(),
              &PeakDef::lessThanByMeanShrdPtr );
  }//while( !candidates.empty() )
  
  const auto detResolution = PeakFitUtils::coarse_resolution_from_peaks( fitpeakvec );
  if( detResolution == PeakFitUtils::CoarseResolutionType::High )
    fitpeakvec = filter_anomolous_width_peaks_highres( meas, fitpeakvec );
  
  return fitpeakvec;
}//search_for_peaks_singlethread(...)

  
vector<std::shared_ptr<const PeakDef> > search_for_peaks(
                              const std::shared_ptr<const Measurement> meas,
                              const std::shared_ptr<const DetectorPeakResponse> drf,
                              std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > origpeaks,
                              const bool singleThreaded,
                              const bool isHPGe )
{
  vector<std::shared_ptr<const PeakDef> > answer;
  
  if( singleThreaded )
    answer = search_for_peaks_singlethread( meas, drf, origpeaks, isHPGe );
  else
    answer = search_for_peaks_multithread( meas, drf, origpeaks, isHPGe );
  
  return answer;
}
  
}//namespace ExperimentalAutomatedPeakSearch


SavitzyGolayCoeffs::SavitzyGolayCoeffs( int bins_left,
                                       int bins_right,
                                       int order,
                                       int derivative )
: num_left( bins_left ),
num_right( bins_right ),
polynomial_order( order ),
ld( derivative ),
coeffs( bins_left + bins_right + 1, 0.0 )
{
  //Savitzy-Golay filter: smoothes data while preserving up to the m'th moment
  //20100315: implemented loosely based on section 14.9 of Numerical Recipes
  //Coefficients stored in order: [lnl, -nl+1, ..., 0, 1, ..., nr]
  //  20111203 - Function adapted from SrbReportsLib
  
  using namespace boost::numeric::ublas;
  
  if( num_left<0 || num_right<0
     || ld>polynomial_order
     || (num_left+num_right)<polynomial_order )
    throw runtime_error( "SavitzyGolayCoeffs(...)\n\tInvalid Input" );
  
  matrix<double> a(polynomial_order+1, polynomial_order+1);
  
  std::vector<float> b(polynomial_order+1, 0.0);
  
  for( int ipj=0; ipj <= (polynomial_order<<1); ++ipj )
  {
    double sum = ( ipj ? 0.0 : 1.0 );
    for( int k=1; k<=num_right; ++k )
      sum += pow( double(k), double(ipj) );
    for( int k=1; k<=num_left; ++k )
      sum += pow( double(-k), double(ipj) );
    const int mm = min(ipj, 2*polynomial_order-ipj);
    for( int imj = -mm; imj<=mm; imj+=2 )
      a((ipj+imj)/2,(ipj-imj)/2) = sum;
  }//for( loop over ipj )
  
  permutation_matrix<std::size_t> pm( a.size1() );
  const size_t res = lu_factorize(a,pm);
  if( res != 0 )
  {
    cerr << "SavitzyGolayCoeffs(...)\n\tFailed to invert Matrix" << endl;
    throw std::runtime_error( "Failed to invert Matrix" );
  }//if( res != 0 )
  
  matrix<double> inverse( a.size1(), a.size1() );
  inverse.assign( identity_matrix<double>(a.size1()) );
  lu_substitute( a, pm, inverse );
  
  //we need only the n^th row of the inverse matrix
  //  meaning this function is computationally inefficient
  for( int i = 0; i <= polynomial_order; ++i )
    b[i] = inverse(ld,i);
  
  for( int k=-num_left; k<=num_right; ++k )
  {
    double sum = b[0];
    double fac = 1.0;
    for( int mm=1; mm <= polynomial_order; ++mm )
      sum += b[mm]*(fac *= k);
    
    coeffs[k+num_left] = sum;
  }//for( loop over kk )
}//SavitzyGolayCoeffs constructor


void SavitzyGolayCoeffs::smooth( const vector<float> &input,
                                vector<float> &output ) const
{
  smooth( &(input[0]), static_cast<int>(input.size()), output );
}


void SavitzyGolayCoeffs::smooth( const float *input,
                                const int nSamples,
                                vector<float> &output ) const
{
  //Performs the Savitzy-Golay filtering with provided coeffs.
  //This function assumes data is flat on either end of the input
  output.clear();
  output.resize( nSamples );
  const int nCoeffs = static_cast<int>( coeffs.size() ); // Of size `bins_left + bins_right + 1`
  
  if( nSamples < nCoeffs || nSamples==0 )
    throw runtime_error( "SavitzyGolayCoeffs::smooth(...)\n\tInvalid input size" );
  
  for( int pos = 0; pos < nSamples; ++pos )
  {
    double sum = 0.0;
    for( int coeff = 0; coeff < nCoeffs; ++coeff )
    {
      int dataInd = pos - num_left + coeff;
      if( dataInd < 0 )
        dataInd = 0;
      else if( dataInd >= nSamples )
        dataInd = nSamples-1;
      
      sum += (coeffs[coeff] * input[dataInd]);
    }//for( loop over coefficients )
    
    output[pos] = sum;
  }//for( loop over input points )
}//SavitzyGolayCoeffs::smooth(...)

std::vector< std::vector<std::shared_ptr<const PeakDef> > >
causilyDisconnectedPeaks(  const double ncausality,
                         const bool useRoiAsWell,
                         std::vector< std::shared_ptr<const PeakDef> > peaks )
{
  typedef std::shared_ptr<const PeakDef> PeakPtr;
  
  std::vector< std::vector<PeakPtr> > answer;
  if( peaks.empty() )
    return answer;
  
  std::sort( peaks.begin(), peaks.end(), &PeakDef::lessThanByMeanShrdPtr );
  
  answer.push_back( vector<PeakPtr>(1,peaks[0]) );
  
  for( size_t i = 1; i < peaks.size(); ++i )
  {
    const PeakPtr &this_peak = peaks[i];
    
    vector<size_t> subPeaksIndicesAddedTo;
    
    for( size_t j = 0; j < answer.size(); ++j )
    {
      vector<PeakPtr> &subpeaks = answer[j];
      for( size_t k = 0; k < subpeaks.size(); ++k )
      {
        const bool isdiscon = PeakDef::causilyDisconnected( *subpeaks[k],
                                         *this_peak, ncausality, useRoiAsWell );
        if( !isdiscon )
        {
          if( subPeaksIndicesAddedTo.empty() )
            subpeaks.push_back( this_peak );
          subPeaksIndicesAddedTo.push_back( j );
          break;
        }
      }//for( size_t k = 0; k < subpeaks.size(); ++k )
    }//for( size_t j = 0; !connected && j < answer.size(); ++j )
    
    if( subPeaksIndicesAddedTo.empty() )
    {
      answer.push_back( vector<PeakPtr>(1,this_peak) );
    }else if( subPeaksIndicesAddedTo.size() > 1 )
    {
      std::sort( subPeaksIndicesAddedTo.begin(), subPeaksIndicesAddedTo.end() );
      
      vector<PeakPtr> &newcombo = answer[subPeaksIndicesAddedTo[0]];
      
      for( size_t j = 1; j < subPeaksIndicesAddedTo.size(); ++j )
      {
        vector<PeakPtr> &oldpeakvec = answer[subPeaksIndicesAddedTo[j]];
        newcombo.insert( newcombo.end(), oldpeakvec.begin(), oldpeakvec.end() );
      }
      
      std::sort( newcombo.begin(), newcombo.end(), &PeakDef::lessThanByMeanShrdPtr );
      
      for( size_t j = subPeaksIndicesAddedTo.size()-1; j > 0; --j )
        answer.erase( answer.begin() + subPeaksIndicesAddedTo[j] );
    }//if( subPeaksIndicesAddedTo.empty() ) / else
  }//for( size_t i = 1; i < peaks.size(); ++i )
  
  return answer;
}//causilyDisconnectedPeaks(...)

      
void unique_copy_continuum( std::vector<PeakDef> &input_peaks )
{
  map<std::shared_ptr<PeakContinuum>,vector<PeakDef>> contToPeaks;
  for( auto &p : input_peaks )
    contToPeaks[p.continuum()].push_back( p );
        
  for( auto &pp : contToPeaks )
  {
    pp.second[0].makeUniqueNewContinuum();
    auto newcont = pp.second[0].continuum();
    for( size_t i = 1; i < pp.second.size(); ++i )
      pp.second[i].setContinuum( newcont );
  }
        
  input_peaks.clear();
  for( auto &pp : contToPeaks )
  {
    for( auto p : pp.second )
      input_peaks.push_back( p );
  }
  std::sort( begin(input_peaks), end(input_peaks), &PeakDef::lessThanByMean );
}//unique_copy_continuum(...)


std::vector< std::vector<PeakDef> > causilyDisconnectedPeaks( const double x0,
                                                             const double x1,
                                                             const double ncausality,
                                                             const bool useRoiAsWell,
                                                             const std::vector<PeakDef> &input_peaks )
{
  const vector<PeakDef> peaks = peaksInRange( x0, x1, ncausality, input_peaks );
  
  std::vector< std::shared_ptr<const PeakDef> > sharedpeaks( peaks.size() );
  for( size_t i = 0; i < peaks.size(); ++i )
    sharedpeaks[i] = std::make_shared<PeakDef>(peaks[i]);
  
  vector< vector<std::shared_ptr<const PeakDef> > > sharedanswer =
             causilyDisconnectedPeaks(  ncausality, useRoiAsWell, sharedpeaks );

  vector< vector<PeakDef> > answer( sharedanswer.size() );
  for( size_t i = 0; i < sharedanswer.size(); ++i )
  {
    for( size_t j = 0; j < sharedanswer[i].size(); ++j )
      answer[i].push_back( *sharedanswer[i][j] );
  }
  
  return answer;
}//causilyDisconnectedPeaks(...)




void findPeaksInUserRange( double x0, double x1, int nPeaks,
                          MultiPeakInitialGuessMethod method,
                          std::shared_ptr<const Measurement> dataH,
                          std::shared_ptr<const DetectorPeakResponse> detector,
                          const bool isHPGe,
                          vector<std::shared_ptr<PeakDef> > &answer,
                          double &chi2 )
{
  //Current (main) problems with the results of this function
  //  --Results can be catastrophically wrong, depending on exact input range
  
  if( method != FromInputPeaks )
    answer.clear();
  chi2 = std::numeric_limits<double>::max();
  
  if( !dataH || nPeaks<=0 )
    return;
  
  if( x1 < x0 )
    std::swap( x0, x1 );
  
  //Lets estimate initial peak parameters
  const size_t start_channel      = dataH->find_gamma_channel( x0 );
  const size_t end_channel        = dataH->find_gamma_channel( x1 );
  const double areaarea    = dataH->gamma_channels_sum( start_channel, end_channel );
  const double start_range = dataH->gamma_channel_lower( start_channel );
  const double end_range   = dataH->gamma_channel_upper( end_channel ) - DBL_EPSILON;
  
  bool intputSharesContinuum = (method==FromInputPeaks);
  for( size_t i = 0; i < answer.size(); ++i )
    intputSharesContinuum &= (answer[i]->continuum()==answer[0]->continuum());
  
  double p0, p1;
  const size_t nSideBinsToAverage = 3;
  PeakContinuum::OffsetType offsetType = PeakContinuum::Linear;
  if( intputSharesContinuum )
    offsetType = answer[0]->continuum()->type();
  else if( nPeaks > 3 && (end_channel - start_channel) > 20 )  //20 is a WAG
    offsetType = PeakContinuum::Quadratic;
  
  PeakContinuum::eqn_from_offsets( start_channel, end_channel, start_range,
                                  dataH, nSideBinsToAverage, nSideBinsToAverage, p1, p0 );
  
  const PeakDef::SkewType skewType = PeakDef::SkewType::NoSkew;
  
  MultiPeakFitChi2Fcn chi2fcn( nPeaks, dataH, offsetType, skewType, start_channel, end_channel );
  
  //chi2fcn.set_reldiff_punish_start( 2.35482 );
  
  vector<PeakDef> inpeaks;
  double conteqn[6] = { p0, p1, 0.0, 0.0, 0.0, 0.0 };
  double initial_cont_area = PeakContinuum::offset_eqn_integral( conteqn, offsetType, start_range, end_range, start_range );
  
  const double totalpeakarea = (areaarea > initial_cont_area)
  ? (areaarea - initial_cont_area) : areaarea;
  
  switch( method )
  {
    case UniformInitialGuess:
    {
      inpeaks.clear();
      for( int i = 0; i < nPeaks; ++i )
      {
        PeakDef peak;
        
        const double amp = totalpeakarea / nPeaks;
        double mean = x0 + (i+0.5)*(x1-x0)/nPeaks;
        double sigma = 0.25*fabs(x1-x0) / nPeaks;
        
        if( detector && detector->hasResolutionInfo() )
          sigma = detector->peakResolutionSigma( mean );
        
        peak.setAmplitude( amp );
        peak.setMean( mean );
        peak.setSigma( sigma );
        
        inpeaks.push_back( peak );
      }//for( int i = 0; i < nPeaks; ++i )
      
      //First well get a little bit better of a guess, by fixing the continuum
      //  to the estimate that it should just touch on each end of the user
      /// defined ROI.  For single this results in fitting for much to wide of a
      
      ROOT::Minuit2::MnUserParameters inputPrams;
      if( offsetType >= PeakContinuum::Constant )
        //        inputPrams.Add( "P0",  p0, std::max( fabs(0.1*p0), 10.0 ) );
        inputPrams.Add( "P0",  p0 );
      if( offsetType >= PeakContinuum::Linear )
        //        inputPrams.Add( "P1",  p1, std::max( fabs(0.1*p1), 10.0 ) );
        inputPrams.Add( "P1",  p1 );
      
      if( offsetType >= PeakContinuum::Quadratic )
        inputPrams.Add( "P2",  0.0 );
      if( offsetType >= PeakContinuum::Cubic )
        inputPrams.Add( "P3",  0.0 );
      
      for( size_t i = 0; i < inpeaks.size(); ++i )
      {
        const string parnum = std::to_string(i);
        
        if( detector && detector->hasResolutionInfo() )
        {
          const double sigma = inpeaks[i].sigma();
          if( i == 0 )
          {
            inputPrams.Add( "sigma" + parnum, sigma, 0.1*sigma, 0.5*sigma, 2.0*sigma );
          }else if( i == 1 )
          {
            const double lowsigma = detector->peakResolutionSigma( start_range );
            const double highsigma = detector->peakResolutionSigma( end_range );
            //sigma = x[m_numOffset] + ((centroid-m_rangeLow)/(m_highRange-m_rangeLow))*x[m_numOffset+3];
            
            const double sigma_diff = (highsigma - lowsigma);
            cout << "sigma_diff=" << sigma_diff << ", highsigma=" << highsigma << ", lowsigma=" << lowsigma << ", sigma=" << sigma << endl;
            if( highsigma > lowsigma && sigma_diff < sigma )
            {
              inputPrams.Add( "sigma" + parnum, sigma_diff, 0.2*sigma_diff, -0.25*sigma, 4*sigma_diff );
            }else
            {
              inputPrams.Add( "sigma" + parnum, 0.0, 0.1, -0.25*sigma, 0.5*sigma );
            }
          }else
          {
            inputPrams.Add( "sigma" + parnum, -1.0 );
          }
        }else
        {
          if( i == 0 )
            inputPrams.Add( "sigma" + parnum, inpeaks[i].sigma(),
                           0.1*inpeaks[i].sigma(), 0, 0.5*(end_range-start_range) );
          else
            inputPrams.Add( "sigma" + parnum, (i==1 ? 1.0 : -1.0) );
        }//if( detector && detector->hasResolutionInfo() ) / else
        
        //        inputPrams.Add( "mean" + parnum, inpeaks[i].mean() );
        const double delta = end_range - start_range;
        inputPrams.Add( "mean" + parnum,
                       inpeaks[i].mean(), 0.1*delta, start_range, end_range );
        
        //        inputPrams.Add( "amplitude" + parnum,
        //                       inpeaks[i].amplitude(),
        //                       0.25*inpeaks[i].amplitude(), 0.0, areaarea );
        //Lets have the chi2 function fit amplitude from the data
        inputPrams.Add( "amplitude" + parnum, -999.9 );
      }//for( size_t i = 0; i < inpeaks.size(); ++i )
      
      //      cerr << "Initial pre-chi2=" << chi2fcn(inputPrams.Params()) << endl;
      
      ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
      
      
      ROOT::Minuit2::MnStrategy strategy( 0 ); //0 low, 1 medium, >=2 high
      ROOT::Minuit2::MnMinimize fitter( chi2fcn, inputParamState, strategy );
      ROOT::Minuit2::FunctionMinimum minimum = fitter( 2500, 1.0 );
      const vector<double> values = fitter.Parameters().Params();
      const vector<double> errors = fitter.Parameters().Errors();
      chi2fcn.parametersToPeaks( inpeaks, &(values[0]), &errors[0] );
      
      if( minimum.IsValid() && (offsetType >= PeakContinuum::Linear) )
      {
        p0 = values[0];
        p1 = values[1];
      }//if( minimum.IsValid() )
      
      //      cerr << "Final pre-chi2=" << chi2fcn(values) << endl;
      break;
    }//case UniformInitialGuess:
      
    case FromDataInitialGuess:
    {
#define USE_QUICK_PEAK_CANDIDATE 1
      
#if( USE_QUICK_PEAK_CANDIDATE )
      intputSharesContinuum = false;
      std::vector< std::tuple<float,float,float> > candidates;  //{mean,sigma,area}
      secondDerivativePeakCanidates( dataH, isHPGe, start_channel, end_channel, candidates );
      std::sort( begin(candidates), end(candidates),
                []( const tuple<float,float,float> &lhs, const tuple<float,float,float> &rhs) -> bool {
                  return std::get<2>(lhs) > std::get<2>(rhs);
      } );
#else
      typedef std::shared_ptr<PeakDef> PeakPtr;
      const vector<PeakPtr> derivative_peaks
        = secondDerivativePeakCanidatesWithROI( dataH, isHPGe, start_channel, end_channel );
      map<double,PeakPtr> candidates;
      for( const PeakPtr &p : derivative_peaks )
        candidates[-p->amplitude()] = p;
#endif
      
    
      double avrg_width = 0.0, avrg_amp = 0.0;
      for( auto i = begin(candidates); int(inpeaks.size()) < nPeaks && i != end(candidates); ++i )
      {
#if( USE_QUICK_PEAK_CANDIDATE )
        if( get<0>(*i) >= x0 && get<0>(*i) <= x1 )
        {
          avrg_width += get<1>(*i);
          avrg_amp += get<2>(*i);
          inpeaks.emplace_back( get<0>(*i), get<1>(*i), get<2>(*i) );
        }//if( i->second.mean() >= x0 && i->second.mean() <= x1 )
#else
        if( i->second->mean() >= x0 && i->second->mean() <= x1 )
        {
          avrg_width += i->second->sigma();
          avrg_amp += i->second->amplitude();
          inpeaks.push_back( *(i->second) );
        }//if( i->second.mean() >= x0 && i->second.mean() <= x1 )
#endif
      }//for(...)
      
      if( inpeaks.size() )
      {
        avrg_width /= inpeaks.size();
        avrg_amp /= inpeaks.size();
      }else
      {
        avrg_width = 0.25*( fabs(x0 - x1) );
        avrg_amp = 0.5*dataH->gamma_channels_sum( start_channel, end_channel );
      }//if( inpeaks ) / else
      
      for( int i = int(inpeaks.size()); i < nPeaks; ++i )
      {
        PeakDef peak;
        const double amp = totalpeakarea / nPeaks;
        const double mean = x0 + (x1-x0)/(1+inpeaks.size());
        double sigma = avrg_width;
        if( detector && detector->hasResolutionInfo() )
          sigma = detector->peakResolutionSigma( mean );
        
        peak.setMean( mean );
        peak.setSigma( sigma );
        peak.setAmplitude( amp );
        
        inpeaks.push_back( peak );
      }//for( int i = int(inpeaks.size()); i < nPeaks; ++i )
      
      break;
    }//case FromDataInitialGuess:
      
    case FromInputPeaks:
    {
      if( answer.size() != static_cast<size_t>(nPeaks) )
        throw runtime_error( "findPeaksInUserRange: invalid input for "
                            "method=FromInputPeaks" );
      for( size_t i = 0; i < answer.size(); ++i )
        inpeaks.push_back( *answer[i] );
      
      answer.clear();
      break;
    }//case FromInputPeaks:
  }//switch( method )
  
  
  ROOT::Minuit2::MnUserParameters inputPrams;
  
  if( intputSharesContinuum )
  {
    const std::vector<double> vals = inpeaks[0].continuum()->parameters();
    for( size_t i = 0; i < vals.size(); ++i )
    {
      const string name = "P" + std::to_string(i);
      const double val = vals[i];
      double delta;
      switch( i )
      {
        case 0:  delta = 10.0; break;
        case 1:  delta = 10.0; break;
        case 2:  delta = 1.0; break;
        default: delta = 0.25; break;
      }//switch( i )
      
      inputPrams.Add( name, val, std::max( fabs(0.1*val), delta ) );
    }//for( size_t i = 0; i < vals.size(); ++i )
  }else
  {
    //calculate po and p1, relative to the first peaks mean
    if( offsetType >= PeakContinuum::Constant )
      inputPrams.Add( "P0",  p0, std::max( fabs(0.1*p0), 10.0 ) );
    if( offsetType >= PeakContinuum::Linear )
      inputPrams.Add( "P1",  p1, std::max( fabs(0.1*p1), 10.0 ) );
    if( offsetType >= PeakContinuum::Quadratic )
      inputPrams.Add( "P2",  0.0, 1.0 );
    if( offsetType >= PeakContinuum::Cubic )
      inputPrams.Add( "P3",  0.0, 0.25 );
  }//if( intputSharesContinuum ) / else
  
  
  for( size_t i = 0; i < inpeaks.size(); ++i )
  {
    const string istr = std::to_string(i);
    double sigma = inpeaks[i].sigma();
    const bool fixSigma = (method==FromInputPeaks && !inpeaks[i].fitFor(PeakDef::Sigma));
    const bool fixMean = (method==FromInputPeaks && !inpeaks[i].fitFor(PeakDef::Mean));
    const bool fixAmp = (method==FromInputPeaks && !inpeaks[i].fitFor(PeakDef::GaussAmplitude));
    
    double maxsigma = fabs(x1-x0)/nPeaks;
    double minsigma = (isHPGe ? 0.0025 : 0.02) * inpeaks[i].mean();
    
    if( !fixSigma )
    {
      float minw, maxw;
      expected_peak_width_limits( inpeaks[i].mean(), isHPGe, dataH, minw, maxw );

      maxsigma = std::max( maxsigma, double(maxw) );
      minsigma = std::min( minsigma, double(minw) );
      sigma = std::max( sigma, minsigma );
      sigma = std::min( sigma, maxsigma );
    }//if( !fixSigma )
    
    if( fixSigma )
      inputPrams.Add( "sigma" + istr, -sigma );
    else if( i == 0 )
      inputPrams.Add( "sigma" + istr, sigma, 0.1*sigma, minsigma, maxsigma );
    else if( i == 1 )
      inputPrams.Add( "sigma" + istr, 0.0, 0.1*sigma, -0.1*sigma, 0.75*sigma ); //could do somethign better here
    else
      inputPrams.Add( "sigma" + istr, 0.0 );
    
    if( fixMean )
      inputPrams.Add( "mean" + istr, inpeaks[i].mean() );
    else
      inputPrams.Add( "mean" + istr, inpeaks[i].mean(), 0.1*fabs(x1-x0), x0 - 0.1, x1 + 0.1 );
    
    if( fixAmp )
    {
      inputPrams.Add( "amplitude" + istr, inpeaks[i].amplitude() );
    }else
    {
      if( inpeaks[i].amplitude() < 25.0 )
      {
        inputPrams.Add( "amplitude" + istr, inpeaks[i].amplitude(),
                       0.25*inpeaks[i].amplitude(), 0.0, areaarea + 1.0 );
      }else
      {
        inputPrams.Add( "amplitude" + istr, (totalpeakarea / nPeaks),
                       (0.25*totalpeakarea / nPeaks), 0.0, totalpeakarea + 1.0 );
      }
    }//if( fixAmp ) / else
  }//for( size_t i = 0; i < inpeaks.size(); ++i )
  
  ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( chi2fcn, inputParamState, strategy );
  
  unsigned int maxFcnCall = 5000;
  const double tolerance = 0.5; //beThorough ? 0.5 : 0.75*nPeaks;
  
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  //cout << "minimum.NFcn()=" << minimum.NFcn() << endl;
  
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
  
  const ROOT::Minuit2::MnUserParameters fitParams = fitter.Parameters();
  const vector<double> values = fitParams.Params();
  const vector<double> errors = fitParams.Errors();
  
  //Turn off punishment for peaks being to close
  chi2fcn.set_reldiff_punish_weight( 0.0 );
  
  chi2 = chi2fcn.DoEval( &(values[0]) );
  const double dof = chi2fcn.dof();
  
  const double chi2Dof = chi2 / dof;
  vector<PeakDef> peaks;
  chi2fcn.parametersToPeaks( peaks, &values[0], &errors[0] );
  
  for( int i = 0; i < nPeaks; ++i )
  {
    auto p = std::make_shared<PeakDef>( peaks[i] );
    p->set_coefficient( chi2Dof, PeakDef::Chi2DOF );
    answer.emplace_back( std::move(p) );
  }
}//void findPeaksInUserRange( double x0, double x1, int nPeaks )





void findPeaksInUserRange_linsubsolve( double x0, double x1, int nPeaks,
                                      MultiPeakInitialGuessMethod method,
                          std::shared_ptr<const Measurement> dataH,
                          std::shared_ptr<const DetectorPeakResponse> detector,
                                      const bool isHPGe,
                          vector<std::shared_ptr<PeakDef> > &answer,
                          double &chi2 )
{
  //Current (main) problems with the results of this function
  //  --Results can be catostropically wrong, depending on exact input range
  
  if( method != FromInputPeaks )
    answer.clear();
  chi2 = std::numeric_limits<double>::max();
  
  if( !dataH || nPeaks<=0 )
    return;
  
  if( x1 < x0 )
    std::swap( x0, x1 );
  
  //Lets estimate initial peak parameters
  const size_t start_channel      = dataH->find_gamma_channel( x0 );
  const size_t end_channel        = dataH->find_gamma_channel( x1 );
  const double areaarea    = dataH->gamma_channels_sum( start_channel, end_channel );
  //const double start_range = dataH->gamma_channel_lower( start_channel );
  //const double end_range   = dataH->gamma_channel_upper( end_channel ) - DBL_EPSILON;
  
  bool intputSharesContinuum = (method==FromInputPeaks);
  for( size_t i = 0; i < answer.size(); ++i )
    intputSharesContinuum &= (answer[i]->continuum()==answer[0]->continuum());
  
  const size_t roiNumChan = (end_channel > start_channel) ? end_channel - start_channel : size_t(0);
  PeakContinuum::OffsetType offsetType = PeakContinuum::Linear;
  if( intputSharesContinuum )
    offsetType = answer[0]->continuum()->type();
  else if( nPeaks > 2 && roiNumChan > 20 && nPeaks < 4 && roiNumChan < 40 )  //20 is a WAG
    offsetType = PeakContinuum::Quadratic;
  else if( nPeaks >= 4 || roiNumChan >= 40 )
    offsetType = PeakContinuum::Cubic;
  
  //MultiPeakFitChi2Fcn chi2fcn( nPeaks, dataH, offsetType, start_channel, end_channel );
  //chi2fcn.set_reldiff_punish_start( 2.35482 );
  
  const PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;
  
  LinearProblemSubSolveChi2Fcn chi2fcn( nPeaks, dataH, offsetType, skew_type, x0, x1 );
  
  float minw_lower, maxw_lower, minw_upper, maxw_upper;
  expected_peak_width_limits( x0, isHPGe, dataH, minw_lower, maxw_lower );
  expected_peak_width_limits( x1, isHPGe, dataH, minw_upper, maxw_upper );
  
  float minsigma = std::min( minw_lower, minw_upper );
  float maxsigma = std::min( maxw_lower, maxw_upper );
  
  if( detector && detector->hasResolutionInfo() )
  {
    minsigma = std::min( minw_lower, detector->peakResolutionSigma(x0) );
    minsigma = std::min( minw_lower, detector->peakResolutionSigma(x1) );
    maxsigma = std::min( maxw_lower, detector->peakResolutionSigma(x0) );
    maxsigma = std::min( maxw_lower, detector->peakResolutionSigma(x1) );
  }
  
  
  /*
   Want to fit for a number of peaks in the given range.
   All peaks should share a common polynomial continuum.
   Widths of the peaks are tied together
   
   parameters:
   mean0     {mean of first peak}
   mean1     {mean of second peak}
   ...
   meanN     {mean of N'th peak}
   width0    {if npeaks==1, width of that peak; else width at m_lowerROI}
   widthFcn  {multiple of width0, plus fraction of total range time widthFcn}
   Skew0     {depends on skew type
   Skew1
   ...
   */

  vector<PeakDef> inpeaks;
  
  switch( method )
  {
    case UniformInitialGuess:
    {
      inpeaks.clear();
      
      for( int i = 0; i < nPeaks; ++i )
      {
        PeakDef peak;
        
        const double amp = areaarea / nPeaks;
        double mean = x0 + (i+0.5)*(x1-x0)/nPeaks;
        double sigma = 0.25*fabs(x1-x0) / nPeaks;
      
        if( detector && detector->hasResolutionInfo() )
          sigma = detector->peakResolutionSigma( mean );
        
        peak.setAmplitude( amp );
        peak.setMean( mean );
        peak.setSigma( sigma );
        
        inpeaks.push_back( peak );
        //inpeaks.emplace_back( std::move(peak) );
      }//for( int i = 0; i < nPeaks; ++i )
      
      break;
    }//case UniformInitialGuess:
      
    case FromDataInitialGuess:
    {
#if( USE_QUICK_PEAK_CANDIDATE )
      intputSharesContinuum = false;
      std::vector< std::tuple<float,float,float> > candidates;  //{mean,sigma,area}
      secondDerivativePeakCanidates( dataH, isHPGe, start_channel, end_channel, candidates );
      std::sort( begin(candidates), end(candidates),
                []( const tuple<float,float,float> &lhs, const tuple<float,float,float> &rhs) -> bool {
                  return std::get<2>(lhs) > std::get<2>(rhs);
                } );
#else
      typedef std::shared_ptr<PeakDef> PeakPtr;
      const vector<PeakPtr> derivative_peaks
      = secondDerivativePeakCanidatesWithROI( dataH, isHPGe, start_channel, end_channel );
      map<double,PeakPtr> candidates;
      for( const PeakPtr &p : derivative_peaks )
        candidates[-p->amplitude()] = p;
#endif
      
    
      double avrg_width = 0.0, avrg_amp = 0.0;
      
      for( auto i = begin(candidates); int(inpeaks.size()) < nPeaks && i != end(candidates); ++i )
      {
#if( USE_QUICK_PEAK_CANDIDATE )
        if( get<0>(*i) >= x0 && get<0>(*i) <= x1 )
        {
          avrg_width += get<1>(*i);
          avrg_amp += get<2>(*i);
          inpeaks.emplace_back( get<0>(*i), get<1>(*i), get<2>(*i) );
        }//if( i->second.mean() >= x0 && i->second.mean() <= x1 )
#else
        if( i->second->mean() >= x0 && i->second->mean() <= x1 )
        {
          avrg_width += i->second->sigma();
          avrg_amp += i->second->amplitude();
          inpeaks.push_back( *(i->second) );
        }//if( i->second.mean() >= x0 && i->second.mean() <= x1 )
#endif
      }//for(...)
      
      if( inpeaks.size() )
      {
        avrg_width /= inpeaks.size();
        avrg_amp /= inpeaks.size();
      }else
      {
        avrg_width = 0.25*( fabs(x0 - x1) );
        avrg_amp = 0.5*dataH->gamma_channels_sum( start_channel, end_channel );
      }//if( inpeaks ) / else
      
      for( int i = int(inpeaks.size()); i < nPeaks; ++i )
      {
        PeakDef peak;
        const double amp = areaarea / nPeaks;
        const double mean = x0 + (x1-x0)/(1+inpeaks.size());
        double sigma = avrg_width;
        if( detector && detector->hasResolutionInfo() )
          sigma = detector->peakResolutionSigma( mean );
        
        peak.setMean( mean );
        peak.setSigma( sigma );
        peak.setAmplitude( amp );
        
        inpeaks.push_back( peak );
      }//for( int i = int(inpeaks.size()); i < nPeaks; ++i )
      
      break;
    }//case FromDataInitialGuess:
      
    case FromInputPeaks:
    {
      if( answer.size() != static_cast<size_t>(nPeaks) )
        throw runtime_error( "findPeaksInUserRange: invalid input for "
                            "method=FromInputPeaks" );
      for( size_t i = 0; i < answer.size(); ++i )
        inpeaks.push_back( *answer[i] );
      
      answer.clear();
      break;
    }//case FromInputPeaks:
  }//switch( method )
  
  assert( !inpeaks.empty() );
  
  /*
   Want to fit for a number of peaks in the given range.
   All peaks should share a common polynomial continuum.
   Widths of the peaks are tied together
   
   parameters:
   mean0     {mean of first peak}
   mean1     {mean of second peak}
   ...
   meanN     {mean of N'th peak}
   width0    {if npeaks==1, width of that peak; else width at m_lowerROI}
   widthFcn  {multiple of width0, plus fraction of total range time widthFcn}
   */
  
  
  
  ROOT::Minuit2::MnUserParameters inputPrams;
  for( size_t i = 0; i < inpeaks.size(); ++i )
  {
    char name[64];
    snprintf( name, sizeof(name), "Mean%i", static_cast<int>(i) );
    const double mean = inpeaks[i].mean();
    
    if( !inpeaks[i].fitFor(PeakDef::Mean) )
      inputPrams.Add( name, mean );
    else
      inputPrams.Add( name, mean, 0.1*(x1-x0), x0, x1 );
  }//for( const PeakDefShrdPtr &peak : inpeaks )
  
  const float sigma0 = static_cast<float>( inpeaks[0].sigma() );
  inputPrams.Add( "Sigma0", std::min(std::max(sigma0,minsigma), maxsigma), 0.2*(maxsigma - minsigma), minsigma, maxsigma );
  
  if( nPeaks > 1 )
    inputPrams.Add( "SigmaFcn", 0.0, 0.001, -0.10, 0.10 );
  
  LinearProblemSubSolveChi2Fcn::addSkewParameters( inputPrams, skew_type );
  
  ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
  ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
  ROOT::Minuit2::MnMinimize fitter( chi2fcn, inputParamState, strategy );
  
  unsigned int maxFcnCall = 5000;
  const double tolerance = 0.5; //beThorough ? 0.5 : 0.75*nPeaks;
  
  ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
  //cout << "minimum.NFcn()=" << minimum.NFcn() << endl;
  
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
  if( !minimum.IsValid() )
    minimum = fitter( maxFcnCall, tolerance );
  
  const ROOT::Minuit2::MnUserParameters fitParams = fitter.Parameters();
  const vector<double> values = fitParams.Params();
  const vector<double> errors = fitParams.Errors();
  
  //Turn off punishment for peaks being to close
  //chi2fcn.set_reldiff_punish_weight( 0.0 );
  
  chi2 = chi2fcn.DoEval( &(values[0]) );
  const double dof = chi2fcn.dof();
  
  const double chi2Dof = chi2 / dof;
  vector<PeakDef> peaks;
  chi2fcn.parametersToPeaks( peaks, &values[0], &errors[0] );
  
  for( int i = 0; i < nPeaks; ++i )
  {
    auto p = std::make_shared<PeakDef>( peaks[i] );
    p->set_coefficient( chi2Dof, PeakDef::Chi2DOF );
    answer.emplace_back( std::move(p) );
  }
}//void findPeaksInUserRange( double x0, double x1, int nPeaks )





void smoothSpectrum( const std::vector<float> &spectrum, const int side_bins,
                    const int order, const int derivative,
                    std::vector<float> &results )
{
  results.clear();
  SavitzyGolayCoeffs sgcoeffs( side_bins, side_bins, order, derivative );
  sgcoeffs.smooth( spectrum, results );
  
  //  if( derivative )
  //  {
  //    for( int bin = 1; bin <= nbin; ++bin )
  //      results[bin-1] /= pow( dataH->GetBinWidth(bin), static_cast<double>(derivative) );
  //  }//if( derivative )
}//void smoothSpectrum(...)

void smoothSpectrum( std::shared_ptr<const Measurement> dataH,
                    const int side,
                    const int order,
                    const int derivative,
                    vector<float> &results )
{
  results.clear();
  if( !dataH || !dataH->gamma_counts() )
    return;
  
  smoothSpectrum( *dataH->gamma_counts(), side, order, derivative, results );
}//void smoothSpectrum(...)



std::vector<PeakDef> fitPeaksInRange( const double x0,
                                     const double x1,
                                     const double ncausality,
                                     const double stat_threshold,
                                     const double hypothesis_threshold,
                                     std::vector<PeakDef> input_peaks,
                                     std::shared_ptr<const Measurement> data,
                                     const std::vector<PeakDef> &fixedpeaks,
                                     bool amplitudeOnly,
                                     const bool isHPGe )
{
  //20120309: For the Ba133 example spectrum with default settings on my newer
  //          mac book pro, this function takes:
  //Single Thread:                   0.069735s wall, 0.070000s user + 0.000000s system = 0.070000s CPU (100.4%)
  //Multithreaded (phys cores only): 0.022624s wall, 0.080000s user + 0.000000s system = 0.080000s CPU (353.6%)
  //Multithreaded (logical cores)  : 0.019906s wall, 0.100000s user + 0.000000s system = 0.100000s CPU (502.3%)
  
  
  typedef vector<PeakDef> PeakVec;
  typedef PeakVec::iterator PeakVecIter;
  
  if( !data || (x1<x0) )
    return input_peaks;
  
  vector<PeakDef> all_peaks = input_peaks;
  
  if( fixedpeaks.size() )
  {
    all_peaks.insert( all_peaks.end(), fixedpeaks.begin(), fixedpeaks.end() );
    std::sort( all_peaks.begin(), all_peaks.end(), &PeakDef::lessThanByMean );
  }//if( fixedpeaks.size() )
  
  
  vector< PeakVec > seperated_peaks
  = causilyDisconnectedPeaks( x0, x1, ncausality, false, all_peaks );
  vector< PeakVec > fixed_seperated_peaks( seperated_peaks.size() );
  
  //Now go through and remove all the fixed peaks from seperated_peaks and
  //  instead place them in fixed_seperated_peaks
  for( size_t index = 0; index < seperated_peaks.size();   )
  {
    PeakVec &peaks = seperated_peaks[index];
    for( const PeakDef &peak : fixedpeaks )
    {
      PeakVecIter pos = lower_bound( peaks.begin(), peaks.end(),
                                    peak, &PeakDef::lessThanByMean );
      if( (pos!=peaks.end()) && ((*pos)==peak) )
      {
        fixed_seperated_peaks[index].push_back( peak );
        peaks.erase( pos );
      }//if( pos != peaks.end() )
    }//for( const PeakDef &peak : fixedpeaks )
    
    if( peaks.empty() )
    {
      seperated_peaks.erase( seperated_peaks.begin() + index );
      fixed_seperated_peaks.erase( fixed_seperated_peaks.begin() + index );
    }else
    {
      ++index;
    }
  }//for( PeakVec peakvec : seperated_peaks )
  
  //remove the peaks we are fitting for from 'input_peaks'
  //  --we will put them back on at the end and return all peaks, not just
  //    in the X range we are fitting in
  for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  {
    const PeakVec &near_peaks = seperated_peaks[peakn];
    
    for( const PeakDef &peak : near_peaks )
    {
      PeakVecIter pos = lower_bound( input_peaks.begin(), input_peaks.end(),
                                    peak, &PeakDef::lessThanByMean );
      if( (pos!=input_peaks.end()) && ((*pos)==peak) )
        input_peaks.erase( pos );
      else cerr << "Logic error in fitPeaksInRange(...)\n";
    }//for( const PeakDef &peak : near_peaks )
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  
  //  boost::timer::cpu_timer timer;
  
  //Fit each of the ranges
  vector< PeakVec > fit_peak_ranges( seperated_peaks.size() );
  SpecUtilsAsync::ThreadPool threadpool;
  //  vector< boost::function<void()> > fit_jobs( seperated_peaks.size() );
  for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  {
    //    fit_jobs[peakn] =
    threadpool.post( boost::bind( &fitPeaks,
                                 boost::cref(seperated_peaks[peakn]),
                                 stat_threshold,
                                 hypothesis_threshold,
                                 data,
                                 boost::ref( fit_peak_ranges[peakn] ),
                                 boost::cref( fixed_seperated_peaks[peakn] ),
                                 amplitudeOnly,
                                 isHPGe ) );
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  
  threadpool.join();
  //  const bool phys_cores_only = false;
  //  SpecUtils::do_asyncronous_work( fit_jobs, phys_cores_only );
  
  
  //put the fit peaks back into 'input_peaks' so we can return all the peaks
  //  passed in, not just ones in X range of interest
  for( size_t peakn = 0; peakn < fit_peak_ranges.size(); ++peakn )
  {
    const PeakVec &fit_peaks = fit_peak_ranges[peakn];
    input_peaks.insert( input_peaks.end(), fit_peaks.begin(), fit_peaks.end() );
  }//for( size_t peakn = 0; peakn < fit_peak_ranges.size(); ++peakn )
  
  std::sort( input_peaks.begin(), input_peaks.end(), &PeakDef::lessThanByMean );
  
  //Now make sure peaks from two previously causally disconnected regions
  //  didn't migrate towards each other, causing the regions to become
  //  causally connected now
  bool migration = false;
  for( size_t peakn = 1; peakn < fit_peak_ranges.size(); ++peakn )
  {
    if( fit_peak_ranges[peakn-1].empty() || fit_peak_ranges[peakn].empty() )
      continue;
    
    const PeakDef &last_peak = fit_peak_ranges[peakn-1].back();
    const PeakDef &this_peak = fit_peak_ranges[peakn][0];
    migration |= PeakDef::causilyConnected( last_peak, this_peak, ncausality, false );
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  
  if( migration )
  {
#ifndef NDEBUG
    cerr << "fitPeaksInRange(...)\n\tWarning: Migration happened!" << endl;
#endif
    
    return fitPeaksInRange( x0, x1, ncausality,
                           stat_threshold, hypothesis_threshold,
                           input_peaks, data, fixedpeaks, amplitudeOnly, isHPGe );
  }//if( migration )
  
  //  cout << "Fit took: " << timer.format() << endl;
  
  return input_peaks;
}//std::vector<PeakDef> fitPeaksInRange(...)


std::shared_ptr<Measurement> estimateContinuum( std::shared_ptr<const Measurement> data )
{
  if( !data )
    throw runtime_error( "estimateContinuum: invalid data" );
  
  int smoothWindow = kBackSmoothing3;  //can be {3, 5, 7, 9, 11, 13, 15}
  int filterOrder = 6; //can be {2, 4, 6, 8}
  const int numIteration = 125; //can be from 1 to 500, roughle
  const bool compton     = false;
  const bool smoothing   = false;
  const int direction    = kBackIncreasingWindow; //kBackDecreasingWindow
  
  const size_t nchannels = data->num_gamma_channels();
  auto source = std::make_shared<vector<float>>( nchannels ) ;
  
  for( size_t i = 0; i < nchannels; ++i )
    (*source)[i] = data->gamma_channel_content(i);

  calculateContinuum( &(source->operator[](0)), static_cast<int>(nchannels), numIteration, direction,
                     filterOrder, smoothing, smoothWindow, compton );
  auto background = make_shared<Measurement>();
  *background = *data;
  
  assert( background->num_gamma_channels() == nchannels );
  
  background->set_gamma_counts( source, data->live_time(), data->real_time() );
  
  return background;
}//std::shared_ptr<Measurement> estimateContinuum( std::shared_ptr<const Measurement> data )



//chi2_for_region(...): gives the chi2 or a region of data, given
//  the input peaks
double chi2_for_region( const PeakShrdVec &peaks,
                       const std::shared_ptr<const Measurement> &data,
                       const int xlowbin,
                       const int xhighbin )
{
  //xxx - need to implement a ROI for multipeaks
  typedef std::shared_ptr<const PeakDef> PeakShrdPtr;
  if( !data || !data->channel_energies() || data->channel_energies()->empty() )
    throw runtime_error( "chi2_for_region: invalid data" );
  
  double chi2 = 0.0;
  
  typedef map< std::shared_ptr<const PeakContinuum>, PeakShrdVec > ContToPeakMap_t;
  ContToPeakMap_t contToPeakMap;
  
  for( size_t peakn = 0; peakn < peaks.size(); ++peakn )
  {
    const std::shared_ptr<const PeakDef> &peak = peaks[peakn];
    std::shared_ptr<const PeakContinuum> continuum = peak->continuum();
    if( continuum->type() == PeakContinuum::External )
      continuum.reset();
    contToPeakMap[continuum].push_back( peak );
  }//for( size_t peak = 0; peak < m_npeaks; ++peak )
  
  for( const ContToPeakMap_t::value_type &vt : contToPeakMap )
  {
    const PeakShrdVec &peaks = vt.second;
    const std::shared_ptr<const PeakContinuum> continuum = peaks[0]->continuum();
    
    if( continuum->energyRangeDefined() )
    {
      const size_t xlowbin = data->find_gamma_channel( continuum->lowerEnergy() );
      const size_t xhighbin = data->find_gamma_channel( continuum->upperEnergy() );
      
      const size_t nchannel = (xhighbin >= xlowbin) ? (1 + xhighbin - xlowbin) : size_t(0);
      vector<double> gauss_counts( std::max(nchannel, size_t(1)), 0.0 );
      const vector<float> &energies = *data->channel_energies();
      for( size_t i = 0; i < peaks.size(); ++i )
        peaks[i]->gauss_integral( &(energies[xlowbin]), &(gauss_counts[0]), nchannel );
      
      for( size_t i = 0; i < nchannel; ++i )
      {
        const size_t channel = xlowbin + i;
        assert( (channel+1) < energies.size() );
        const double xbinlow = energies[channel];
        const double xbinup = energies[channel+1];
        const double ndata = data->gamma_channel_content(channel);
        const double ncontinuum = continuum->offset_integral(xbinlow, xbinup, data);
        
        const double npeak = gauss_counts[i];
        const double uncert = ndata > MIN_CHANNEL_UNCERT ? sqrt(ndata) : 1.0;
        const double chi = (ndata-ncontinuum-npeak)/uncert;
        
        chi2 += chi*chi;
      }//for( size_t i = 0; i <= nchannel; ++i )
    }else
    {
      // 20231107: I'm not sure if we should ever get here anymore
      map<size_t,double> predicted;
      
      for( const PeakShrdPtr &peak : peaks )
      {
        size_t lower_channel = data->find_gamma_channel( peak->lowerX() );
        size_t upper_channel = data->find_gamma_channel( peak->upperX() - 0.00001 );
        
        for( size_t channel = lower_channel; channel <= upper_channel; ++channel )
        {
          if( !predicted.count(channel) )
            predicted[channel] = 0.0;
          const double xbinlow = data->gamma_channel_lower(channel);
          const double xbinup = data->gamma_channel_upper(channel);
          const double ncontinuum = continuum->offset_integral(xbinlow, xbinup, data);
          predicted[channel] += ncontinuum + peak->gauss_integral( xbinlow, xbinup );
        }//for( int bin = xlowbin; bin <= xhighbin; ++bin )
      }//for( const PeakShrdPtr &peak : peaks )
      
      for( map<size_t,double>::iterator i = predicted.begin(); i != predicted.end(); ++i )
      {
        const double peakarea = i->second;
        const double ndata = data->gamma_channel_content( i->first );
        const double uncert = ndata > MIN_CHANNEL_UNCERT ? sqrt(ndata) : 1.0;
        const double chi = (ndata-peakarea)/uncert;
        chi2 += chi*chi;
      }//for( int bin = minbin; bin <= maxbin; ++bin )
    }//if( continuum->energyRangeDefined() ) / else
  }//for( const ContToPeakMap_t::value_type &vt : contToPeakMap )
  
  return chi2;
}//chi2_for_region(...)



PeakShrdVec refitPeaksThatShareROI( const std::shared_ptr<const Measurement> &dataH,
                                   const DetctorPtr &detector,
                                   const PeakShrdVec &inpeaks,
                                   const double meanSigmaVary )
{
  typedef std::shared_ptr<const PeakDef> PeakPtr;
  
  PeakShrdVec answer;
  
  try
  {
    if( inpeaks.empty() )
      return answer;
    
    std::shared_ptr<const PeakContinuum> origCont = inpeaks[0]->continuum();
    
    for( const PeakPtr &p : inpeaks )
      if( origCont != p->continuum() )
        throw runtime_error( "refitPeaksThatShareROI: all input peaks must share a ROI" );
    
    const double lx = origCont->lowerEnergy();
    const double ux = origCont->upperEnergy();
    const double range = ux - lx;
    
    
    PeakDef::SkewType skew_type = LinearProblemSubSolveChi2Fcn::skewTypeFromPrevPeaks( inpeaks );
    
    
    LinearProblemSubSolveChi2Fcn chi2Fcn( inpeaks, dataH, origCont->type(), skew_type, lx, ux );
    
    const size_t mid_channel = dataH->find_gamma_channel( lx + 0.5*range );
    double minsigma = dataH->gamma_channel_width( mid_channel );
    if( detector && detector->hasResolutionInfo() )
      minsigma = 0.75*detector->peakResolutionSigma( origCont->lowerEnergy() );
    
    int nFitWidth = 0, nFitEnergy = 0, nSkewFitPars = 0;
    ROOT::Minuit2::MnUserParameters params;
    for( size_t i = 0; i < inpeaks.size(); ++i )
    {
      if( !inpeaks[i]->gausPeak() )
        throw runtime_error( "Somehow a data defined peak made it into "
                            "refitPeaksThatShareROI; please let "
                            "InterSpec@sandia.gov know about this" );
      
      char name[64];
      snprintf( name, sizeof(name), "Mean%i", static_cast<int>(i) );
      const double mean = inpeaks[i]->mean();
      const double sigma = inpeaks[i]->sigma(); //Will throw exception if peak not Gaussian defined - wanted behavior
      minsigma = std::min( minsigma, 0.9*sigma );
      
      if( !inpeaks[i]->fitFor(PeakDef::Mean) || meanSigmaVary==0.0 )
        params.Add( name, mean );
      else if( meanSigmaVary > 0.0 )
        params.Add( name, mean, 0.1*sigma, mean-meanSigmaVary*sigma, mean+meanSigmaVary*sigma );
      else
        params.Add( name, mean, 0.1*sigma, lx, ux );
      
      nFitWidth += inpeaks[i]->fitFor(PeakDef::Sigma);
      nFitEnergy += inpeaks[i]->fitFor(PeakDef::Mean);
      
      for( int skew_par = 0; skew_par < PeakDef::num_skew_parameters(skew_type); ++skew_par )
      {
        const auto ct = PeakDef::CoefficientType(PeakDef::CoefficientType::SkewPar0 + skew_par);
        nSkewFitPars += inpeaks[i]->fitFor(ct);
      }
    }//for( const PeakDefShrdPtr &peak : inpeaks )
    
      
    const double sigma0 = inpeaks[0]->sigma();
    if( nFitWidth == 0  )
      params.Add( "Sigma0", sigma0 );
    else
      params.Add( "Sigma0", sigma0, 0.1*sigma0, minsigma, 0.75*range );
    
    if( nFitWidth > 1 )
      params.Add( "SigmaFcn", 0.0, 0.001, -0.10, 0.10 );
    else if( inpeaks.size() > 1 )
      params.Add( "SigmaFcn", 0.0 );
    
    if( skew_type != PeakDef::SkewType::NoSkew )
    {
      if( meanSigmaVary > 0.3 )
        LinearProblemSubSolveChi2Fcn::addSkewParameters( params, skew_type );
      else
        LinearProblemSubSolveChi2Fcn::addSkewParameters( params, skew_type, inpeaks );
    }//if( skew_type != PeakDef::SkewType::NoSkew )
    
    
    if( (nFitWidth == 0) && (nFitEnergy == 0) && (nSkewFitPars == 0) )
    {
      // Nothing to do here; LinearProblemSubSolveChi2Fcn::parametersToPeaks(...) will do all work
    }else
    {
      ROOT::Minuit2::MnUserParameterState inputParamState( params );
      ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
      
      ROOT::Minuit2::CombinedMinimizer fitter;
      unsigned int maxFcnCall = 0;
      const double tolerance = 0.01;
      
      auto minimum = fitter.Minimize( chi2Fcn, params, strategy, maxFcnCall, tolerance );
      
      params = minimum.UserState().Parameters();
    }//if( no non-linear fit parameters ) / else
    
    const vector<double> pars = params.Params();
    const vector<double> errors = params.Errors();
      
    vector<PeakDef> fitpeaks;
    double chi2 = chi2Fcn.parametersToPeaks( fitpeaks, &pars[0], &errors[0] );
    double chi2Dof = chi2 / chi2Fcn.dof();
    
    if( fitpeaks.size() != inpeaks.size() )
      throw std::logic_error( "refitPeaksThatShareROI: invalid number of result peaks" );
    
    static int ntimesmessages = 0;
    if( ntimesmessages++ < 3 )
    {
      cerr << "The fit yielded chi2Dof=" << chi2Dof << endl;
      cerr << "Need to check here if the fit was any good!!!" << endl;
    }
    
    for( size_t i = 0; i < fitpeaks.size(); ++i )
    {
      fitpeaks[i].set_coefficient( chi2Dof, PeakDef::Chi2DOF );
      answer.push_back( std::make_shared<PeakDef>(fitpeaks[i]) );
    }//for( size_t i = 0; i < fitpeaks.size(); ++i )
    
    //now we need to go through and make sure the peaks we're adding are both
    //  significant, and improve the chi2/dof from before the fit.
    const int lower_channel = static_cast<int>( dataH->find_gamma_channel( lx ) );
    const int upper_channel = static_cast<int>( dataH->find_gamma_channel( ux ) );
    const int nbin = (upper_channel > lower_channel) ? (upper_channel - lower_channel) : 1;
    const double prechi2Dof = chi2_for_region( inpeaks, dataH, lower_channel, upper_channel ) / nbin;
    const double postchi2Dof = chi2_for_region( answer, dataH, lower_channel, upper_channel ) / nbin;
    
    if( ntimesmessages++ < 4 )
    {
      cerr << "prechi2Dof=" << prechi2Dof << endl;
      cerr << "postchi2Dof=" << postchi2Dof << endl;
    }
    
    if( prechi2Dof < postchi2Dof )
    {
      answer.clear();
      
      if( true )
      {
        const double ncausalitysigma = 0.0;
        const double stat_threshold  = 0.0;
        const double hypothesis_threshold = 0.0;
        
        const bool isRefit = true;
        
        vector<PeakDef> input_peaks, fixed_peaks;
        for( const auto &p : inpeaks )
          input_peaks.push_back( *p );
        
        const auto resType = PeakFitUtils::coarse_resolution_from_peaks(inpeaks);
        const bool isHPGe = (resType == PeakFitUtils::CoarseResolutionType::High);
        
        
        vector<PeakDef> output_peak = fitPeaksInRange( lx, ux, ncausalitysigma, stat_threshold,
                                     hypothesis_threshold, input_peaks, dataH,
                                     fixed_peaks, isRefit, isHPGe );
        
        if( output_peak.size() == inpeaks.size() )
        {
          vector<shared_ptr<const PeakDef>> refit_peaks;
          for( const auto &p : output_peak )
            refit_peaks.push_back( make_shared<PeakDef>(p) );
          
          const double refit_chi2Dof = chi2_for_region( refit_peaks, dataH, lower_channel, upper_channel ) / nbin;
          //cout << "refit_chi2Dof=" << refit_chi2Dof << ", prechi2Dof=" << prechi2Dof << ", postchi2Dof=" << postchi2Dof << endl;
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
    
      
    for( const PeakPtr peak : answer )
    {
      const double mean = peak->mean();
      const double sigma = peak->sigma();
      const double data = gamma_integral( dataH, mean-0.5*sigma, mean+0.5*sigma );
      const double area = peak->gauss_integral( mean-0.5*sigma, mean+0.5*sigma );
      
      if( area < 5.0 || area < sqrt(data) )
      {
        answer.clear();
        break;
      }//if( area < 5.0 || area < sqrt(data) )
      
      if( answer.size() == 1 )
      {
        static int ntimes = 0;
        if( ntimes++ < 3 )
          cerr << "refitPeaksThatShareROI: Should perform some additional chi2 checks for the "
          "case answer.first.size() == 1" << endl;
      }//if( answer.first.size() == 1 )
    }//for( const PeakDefShrdPtr peak : answer.first )
    
  }catch( std::exception & )
  {
    answer.clear();
    cerr << "refitPeaksThatShareROI: failed to find a better fit" << endl;
  }
  
  return answer;
}//PeakShrdVec refitPeaksThatShareROI(...)


double evaluate_chi2dof_for_range( const std::vector<PeakDef> &peaks,
                                  const std::shared_ptr<const Measurement> &dataH,
                                  const double startx,
                                  const double endx )
{
  const size_t lowerchannel = dataH->find_gamma_channel( startx );
  const size_t upperchannel = dataH->find_gamma_channel( endx );
  
  double chi2 = 0;
  for( size_t channel = lowerchannel; channel <= upperchannel; ++channel )
  {
    const double x0 = dataH->gamma_channel_lower( channel );
    const double x1 = dataH->gamma_channel_lower( channel + 1 );
    
    const float y = dataH->gamma_channel_content( channel );
    
    double y_pred = 0.0;
    
    set<const PeakContinuum *> continuums;
    
    for( size_t j = 0; j < peaks.size(); ++j )
    {
      const PeakContinuum * const contptr = peaks[j].continuum().get();
      
      if( x1 < peaks[j].lowerX() || x0 > peaks[j].upperX() )
        continue;
      
      y_pred += peaks[j].gauss_integral(x0, x1);
      
      if( !continuums.count( contptr ) )
      {
        continuums.insert( contptr );
        y_pred += contptr->offset_integral(x0, x1, dataH);
      }
    }//for( size_t j = 0; j < peaks.size(); ++j )
    
    if( y_pred < 0.0 )
      y_pred = 0.0;
    
    const double uncert = (y > MIN_CHANNEL_UNCERT ? sqrt(y) : 1.0);
    chi2 += std::pow( (y_pred - y) / uncert, 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2 / (1 + upperchannel - lowerchannel);
}//double evaluate_chi2dof_for_range(...);


//refit_for_new_roi(): resultChi2Dof will be DBL_MAX and resultPeaks empty upon
//  error.
void refit_for_new_roi( std::vector< std::shared_ptr<const PeakDef> > originalPeaks,
                       const std::shared_ptr<const Measurement> &dataH,
                       const double new_lower_roi,
                       const double new_roi_upper,
                       double &resultChi2Dof,
                       std::vector<PeakDef> &resultPeaks )
{
  try
  {
    const PeakContinuum::OffsetType offset = originalPeaks[0]->continuum()->type();
    
    const PeakDef::SkewType skew_type = LinearProblemSubSolveChi2Fcn::skewTypeFromPrevPeaks( originalPeaks );
    
    LinearProblemSubSolveChi2Fcn chi2Fcn( originalPeaks, dataH, offset, skew_type,
                                         new_lower_roi, new_roi_upper );
    
    //If there are a ton of peaks in this ROI, a compleete refit will
    //  be to slow, so only to the complete re-fit for 1 or 2 peaks.
    const size_t npeaks = originalPeaks.size();
    
    ROOT::Minuit2::MnUserParameters params;
    
    for( size_t i = 0; i < originalPeaks.size(); ++i )
    {
      char name[64];
      const std::shared_ptr<const PeakDef> &peak = originalPeaks[i];
      snprintf( name, sizeof(name), "Mean%i", static_cast<int>(i) );
      const double mean = peak->mean();
      const double sigma = peak->sigma();
      
      const double meanmin = std::max( new_lower_roi, mean-1.0*sigma );
      const double meanmax = std::min( new_roi_upper, mean+1.0*sigma );
      
      if( !peak->fitFor(PeakDef::Mean) )
        params.Add( name, mean );
      else
        params.Add( name, mean, 0.1*sigma, meanmin, meanmax );
    }//for( const PeakDefShrdPtr &peak : peaksToReFit )
    
    const double sigma0 = originalPeaks[0]->sigma();
    
    params.Add( "Sigma0", sigma0, 0.1*sigma0, 0.25*sigma0, 2.0*sigma0 );
    if( npeaks > 1 )
      params.Add( "SigmaFcn", 0.0, 0.001, -0.10, 0.10 );
    
    LinearProblemSubSolveChi2Fcn::addSkewParameters( params, skew_type, originalPeaks );
    
    ROOT::Minuit2::MnUserParameterState inputParamState( params );
    const ROOT::Minuit2::MnStrategy strategy( 2 );
    
    ROOT::Minuit2::CombinedMinimizer fitter;
    ROOT::Minuit2::FunctionMinimum minimum
    = fitter.Minimize( chi2Fcn, params, strategy, 0, 0.01 );
    
    params = minimum.UserState().Parameters();
    const vector<double> pars = params.Params();
    const vector<double> errors = params.Errors();
    
    const double newChi2 = chi2Fcn.parametersToPeaks( resultPeaks, &pars[0], &errors[0] );
    
    resultChi2Dof = newChi2 / chi2Fcn.dof();
  }catch( std::exception & )
  {
    resultPeaks.clear();
    resultChi2Dof = DBL_MAX;
  }//try / catch
}//void refit_for_new_roi(...)



//Throws exception on error
void find_roi_for_2nd_deriv_candidate(
                                      double &lowerEnengy,
                                      double &upperEnergy,
                                      const float peakmean,
                                      const std::shared_ptr<const Measurement> &data,
                                      const bool isHPGe )
{
  if( !data || !data->num_gamma_channels() )
    throw runtime_error( "find_roi_for_2nd_deriv_candidate: invalid input" );
  
  const size_t nchannel = data->num_gamma_channels();
  
  const size_t meanchannel = data->find_gamma_channel( peakmean );
  
  //Changing spoothing to be based in middle channel, not the channel of the
  //  candidate mean, to be more consistent with
  //  const float binwidth = data->gamma_channel_width( meanchannel );
  //  int side_bins = isHPGe ? 4 : std::max( 5, static_cast<int>( floor(0.022*peakmean/binwidth+0.5) ) );
  //  const size_t midchannel = data->num_gamma_channels() / 2;
  //  const float midenergy = data->gamma_channel_lower( midchannel );
  //  const float binwidth = data->gamma_channel_width( midchannel );
  //  const int side_bins = isHPGe ? 4 : std::max( 5, static_cast<int>( floor(0.022*peakmean/binwidth+0.5) ) );
  
  //   const int order = isHPGe ? 3 : 2;
  
  //The below should be same as in secondDerivativePeakCanidatesWithROI(...)
  const size_t midbin = data->num_gamma_channels() / 2;// (start_channel + end_channel) / 2;
  const float midenergy = data->gamma_channel_center( midbin );
  const float midbinwidth = data->gamma_channel_width( midbin );
  
  const int order = isHPGe ? 3 : 2;
  const size_t side_bins = isHPGe ? 4 : std::max( size_t(5), static_cast<size_t>( 0.022f*midenergy/midbinwidth + 0.5f ) );
  

  vector<float> smoothed, second_deriv;
  smoothSpectrum( data, static_cast<int>(side_bins), order, 0, smoothed );
  smoothSpectrum( data, static_cast<int>(side_bins), order, 2, second_deriv );
  
  
  //This next part should be kept the same as in
  //  secondDerivativePeakCanidatesWithROI(...) during development.
  if( !isHPGe && side_bins > 5 && nchannel >= 512 )
  {
    //We should also have a minbimum statistics requirment here.
    
    const size_t index = (nchannel/15);
    vector<float> second_deriv_lower, smoothed_lower;
    smoothSpectrum( data, 4, order, 0, smoothed_lower );
    smoothSpectrum( data, 4, order, 2, second_deriv_lower );
    
    for( size_t i = 0; i < (index-side_bins); ++i )
    {
      smoothed[i] = smoothed_lower[i];
      second_deriv[i] = second_deriv_lower[i];
    }
    
    //transition over 'side_bins' between the smoothings.
    for( size_t i = 0; i < side_bins; ++i )
    {
      const float factor = float(i+1) / float(side_bins+1);
      const size_t current = index - side_bins + i;
      second_deriv[current] = factor*second_deriv[current]
                              + (1.0f-factor)*second_deriv_lower[current];
      smoothed[current] = factor*smoothed[current]
                          + (1.0f-factor)*smoothed_lower[current];
    }
    
//    cout << "Transition occurs at " << data->gamma_channel_center(index) << " kev" << endl;
  }//if( !isHPGe )
  
  
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  {
    {
      ofstream secondfile( "secondderivROI.csv" ), smoothfile( "smoothedROI.csv" );
      secondfile << "remark: livetime " << data->live_time()/500.0f << "s" << endl;
      smoothfile << "remark: livetime " << data->live_time() << "s" << endl;
      secondfile << "Energy,Counts" << endl;
      smoothfile << "Energy,Counts" << endl;
      for( size_t bin = 0; bin < nchannel; ++bin )
      {
        secondfile << data->gamma_channel_lower(bin) << "," << second_deriv[bin] << endl;
        smoothfile << data->gamma_channel_lower(bin) << "," << smoothed[bin] << endl;
      }
    }
    DebugLog(cout) << "Made secondderivROI.csv and smoothedROI.csv for mean " << peakmean << "\n"
                   << "side_bins=" << side_bins << ", order=" << order << "\n";
  }
#endif
  
  if( second_deriv[meanchannel] > 0.0 )
    throw runtime_error( "find_roi_for_2nd_deriv_candidate: invalid mean" );
  
  //In order to allow small flucations right around transitions from above or
  //  below y=zero to the other side, we will enforce 'nFluxuateBin' having
  //  to be either above or below zero before we've declare a transisition has
  //  happened.
  //This number should probably be kept the same as in
  //  secondDerivativePeakCanidatesWithROI(...)
  const float meanuncert = 1.0f / sqrt( smoothed[meanchannel] );
  const size_t nFluxuateBin = (meanuncert < 0.05f) ? ((meanuncert < 0.015f) ? 2 : 3) : 4;
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  DebugLog(cout) << "mean=" << peakmean << ", uncert=" << meanuncert << "\n";
#endif
  
  assert( nFluxuateBin >= 1 );
  
  size_t firstnegzero = meanchannel;
  while( firstnegzero > nFluxuateBin )
  {
    bool above = true;
    for( size_t i = 0; i < nFluxuateBin; ++i )
      above &= (second_deriv[firstnegzero-i] >= 0.0f);
    if( above )
      break;
    --firstnegzero;
  }
  
  size_t secondnegzero = firstnegzero;
  while( secondnegzero > nFluxuateBin )
  {
    bool below = true;
    for( size_t i = 0; i < nFluxuateBin; ++i )
      below &= (second_deriv[secondnegzero-i] < 0.0f);
    if( below )
      break;
    --secondnegzero;
  }
  
  size_t firstposzero = meanchannel;
  while( firstposzero < (nchannel-nFluxuateBin) )
  {
    bool above = true;
    for( size_t i = 0; i < nFluxuateBin; ++i )
      above &= (second_deriv[firstposzero+i] >= 0.0f);
    if( above )
      break;
    ++firstposzero;
  }
  
  size_t secondposzero = firstposzero;
  while( secondposzero < (nchannel-nFluxuateBin) )
  {
    bool below = true;
    for( size_t i = 0; i < nFluxuateBin; ++i )
      below &= (second_deriv[secondposzero+i] < 0.0f);
    if( below )
      break;
    ++secondposzero;
  }
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 1 )
  DebugLog(cout) << "First negzero at " << data->gamma_channel_lower(firstnegzero)
  << ", secondnegzero at " << data->gamma_channel_lower(secondnegzero)
  << ", firstposzero at " << data->gamma_channel_lower(firstposzero)
  << ", secondposzero at " << data->gamma_channel_lower(secondposzero)
  << ", meanbin at " << data->gamma_channel_lower(meanchannel)
  <<"\n";
#endif
  
  
  for( size_t lowchannel = secondnegzero; lowchannel < firstnegzero; ++lowchannel)
  {
    for( size_t highchannel = secondposzero; highchannel > firstposzero; --highchannel )
    {
      const float y0 = smoothed[lowchannel];
      const float y1 = smoothed[highchannel];
      
      const float m = (y1 - y0) / (highchannel - lowchannel);
      const float b = y1 - m*highchannel;
      
      bool belowdata = true;
      for( size_t i = secondnegzero; belowdata && i <= secondposzero; ++i )
      {
        if( i != lowchannel && i != highchannel )
          belowdata &= (smoothed[i] >= (i*m+b));
      }
      
      if( belowdata )
      {
        const float lintersection = data->gamma_channel_center( lowchannel );
        const float rintersection = data->gamma_channel_center( highchannel );
        
        lowerEnengy = lintersection;
        upperEnergy = rintersection;
      
        float avrguncert = 0.0f;
        for( size_t i = lowchannel; i <= highchannel; ++i )
          avrguncert += (1.0f / std::sqrt( std::max( smoothed[i], 1.0f ) ));
        avrguncert /= (highchannel - lowchannel + 1);
        
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 1 )
        DebugLog(cout) << "avrguncert=" << avrguncert << " for " << peakmean << "\n";
#endif
      
        //try to extend the ROI
        //  the 0.1, 0.35, and 1.25 are all chosen haphazerdly, with no real
        //  inteligence, and not a ton of testing
        const bool highstats = (avrguncert < 0.1f);
        const float max_roi_extend_frac = highstats ? 0.35f : 1.0f;
        const float max_delta_devaition = 5.0f;
        
        const float llength = peakmean - lintersection;
        const float rlength = upperEnergy - peakmean;
        const float minroienergy = lintersection - max_roi_extend_frac*llength;
        const float maxroienergy = rintersection + max_roi_extend_frac*rlength;
        const size_t minroibin = std::max( size_t(1), data->find_gamma_channel( minroienergy ) );
        const size_t maxroibin = data->find_gamma_channel( maxroienergy );
        
//        cout << "lowchannel=" << lowchannel << ", highchannel=" << highchannel
//        << ", lowerEnengy=" << lowerEnengy << ", upperEnergy=" << upperEnergy << endl;
        
        //If sum of area between continuum and data is negative (line above
        //  data), keep trying to extend the ROI
        //If line is below data on average, keep going until sum of area
        //  between continuum and data, and the sum of data diaviate by
        float diffsum = 0.0f, datasum = 0.0f;
        for( size_t i = lowchannel; i >= minroibin; --i )
        {
          const float dataval = smoothed[i];
          const float lineval = i*m + b;
          datasum += dataval;
          diffsum += (dataval - lineval);
          const float dev = diffsum/sqrt(std::max(datasum,1.0f));

          lowchannel = i + 1;
          lowerEnengy = data->gamma_channel_center( lowchannel );

#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 1 )
          DebugLog(cout) << "Extending lower edge to " << lowerEnengy
               << ", dev=" << dev << "\n";
#endif
          
          if( i != lowchannel && fabs(dev) > max_delta_devaition )
            break;
        }//for( size_t i = lowchannel; i > minroibin; --i )
        
        diffsum = datasum = 0.0f;
        for( size_t i = highchannel; i <= maxroibin; ++i )
        {
          const float dataval = data->gamma_channel_content( i ); //smoothed[i];
          const float lineval = i*m + b;
          datasum += dataval;
          diffsum += (dataval - lineval);
          const float dev = diffsum/sqrt(std::max(datasum,1.0f));

          highchannel = i - 1;
          upperEnergy = data->gamma_channel_center( highchannel );
          
          if( i != highchannel && fabs(dev) > max_delta_devaition )
            break;
        }//for( size_t i = highchannel; i <= maxroibin; ++i )
        
        //Since the detector has a 'turn-on' for low energies (where the
        //  spectrum is ramping up, as efficinecy increases with energy), we
        //  need to make sure we dont try to set the ROI to lower than the part
        //  of the spectrum that meaningfully carries peak information.  We'll
        //  do this by seeing if the ROI lower goes below the spectroscopic
        //  extent (we'l aslo only check this when lowerEnergy is below 100
        //  keV).
        if( lowerEnengy < 100.0 )
        {
          size_t speclower = 0, specupper = 0;
          ExperimentalPeakSearch::find_spectroscopic_extent( data, speclower, specupper );
          if( lowchannel < speclower )
          {
            //The lowere energy portion of the spectrum behaves so badly, we'll
            //  just use where the second derivative crosses zero
            lowchannel = std::max( lowchannel, firstnegzero );
            lowerEnengy = data->gamma_channel_center( lowchannel );
          }
        }//if( lowerEnengy < 100.0 )
        
        return;
      }//if( belowdata )
    }//for( size_t highchannel = secondposzero; highchannel > firstposzero; --highchannel )
  }//for( size_t lowchannel = secondnegzero; lowchannel < firstnegzero; ++lowchannel)
  
  throw runtime_error( "find_roi_for_2nd_deriv_candidate: "
                      "failed to fine line for ROI definition" );
}//void find_roi_for_2nd_deriv_candidate(...)


void expected_peak_width_limits( const float energy,
                                 const bool highres,
                                 const std::shared_ptr<const SpecUtils::Measurement> &meas,
                                 float &min_sigma_width_kev,
                                 float &max_sigma_width_kev )
{
  //For this function I took approximately the best and worst resolution
  //  detector response functions for lowres detectors in GADRAS, and this
  //  function just gives a mutliple of these values.
  //The arrays in this function were generated using the
  //  print_detector_sigma_range() function in developcode.cpp.
  const float max_width_multple = highres ? 4.0 : 3.0f;
  const float min_width_multiple = highres ? (energy < 50.0f ? 0.25f : 0.35f) : 0.5f;
  
  const size_t nenergies = 32;
  const float energies[32]
  = { 10.0f, 25.0f, 35.0f, 45.0f, 55.0f, 70.0f, 90.0f, 120.0f,
    170.0f, 250.0f, 350.0f, 450.0f, 550.0f, 650.0f, 750.0f, 850.0f,
    950.0f, 1100.0f, 1300.0f, 1500.0f, 1700.0f, 1900.0f, 2200.0f, 2600.0f,
    3000.0f, 4000.0f, 5000.0f, 7000.0f, 10000.0f, 15000.0f, 30000.0f, 60000.0f };

  const float lowres_smallest_expected_sigma[32]
  = { 0.827749f, 1.5699f, 1.78745f, 1.98339f, 2.16319f, 2.40991f, 2.70684f, 3.10335f,
    3.67698f, 4.45373f, 5.27675f, 5.99548f, 6.64232f, 6.9945f, 7.19707f, 7.3795f,
    7.5455f, 7.77002f, 8.03401f, 8.26726f, 8.47683f, 8.66751f, 8.92541f, 9.22865f,
    9.49659f, 10.059f, 10.5181f, 11.2503f, 12.0821f, 13.1027f, 15.0511f, 17.2892f };
  const float lowres_largest_expected_sigma[32]
  = { 3.50046f, 4.64015f, 5.4903f, 6.22541f, 6.88245f, 7.76445f, 8.80406f, 10.1661f,
    12.1f, 14.6734f, 17.3618f, 19.6865f, 21.7642f, 23.6602f, 25.4151f, 27.0565f,
    28.6038f, 30.7793f, 33.4606f, 35.9424f, 38.7333f, 41.5448f, 45.5646f, 50.6215f,
    55.3972f, 66.4048f, 76.428f, 94.4742f, 118.277f, 152.7f, 236.313f, 365.71f };

  const float highres_smallest_expected_sigma[32]
  = { 0.42612f, 0.428721f, 0.430567f, 0.432469f, 0.434415f, 0.437396f,
      0.441456f, 0.447676f, 0.45825f, 0.47541f, 0.496911f, 0.518233f, 0.539259f,
      0.559936f, 0.603397f, 0.647207f, 0.688801f, 0.747736f, 0.821064f,
      0.889569f, 0.945747f, 0.975664f, 1.01655f, 1.06523f, 1.10877f, 1.20178f,
      1.27927f, 1.40565f, 1.55328f, 1.74002f, 2.11273f, 2.56526f };
  const float highres_largest_expected_sigma[32]
  = { 0.752363f, 0.752363f, 0.752363f, 0.752363f, 0.752363f, 0.752363f,
      0.752363f, 0.752363f, 0.752363f, 0.752363f, 0.752363f, 0.752363f,
      0.752363f, 0.752363f, 0.772945f, 0.805539f, 0.835656f, 0.885211f,
      0.974007f, 1.05712f, 1.13561f, 1.21024f, 1.31615f, 1.44817f, 1.57175f,
      1.85301f, 2.10539f, 2.55241f, 3.13033f, 3.9478f, 5.86967f, 8.72714f };
  
  
  const float *beginx = energies;
  const float *endx = beginx + nenergies;
  const float *pos = std::lower_bound( beginx, endx, energy );
  
  //I dont expect that energy range will ever be invalid, so not taking to much
  //  care to deal with this possiblity (other than prevent a crash)
  if( pos >= (endx-1) )
    pos = endx - 2;
  
  const size_t index = pos - beginx;
  const float prev_energy = energies[index];
  const float next_energy = energies[index+1];
  const float frac = (energy - prev_energy) / (next_energy - prev_energy);
  
  const float *smallest = highres ? highres_smallest_expected_sigma + 0
                                  : lowres_smallest_expected_sigma + 0;
  const float *largest  = highres ? highres_largest_expected_sigma + 0
                                  : lowres_largest_expected_sigma + 0;
  
  const float deltalarge = largest[index+1] - largest[index];
  const float deltasmall = smallest[index+1] - smallest[index];
  
  const float largerres = largest[index] + frac*deltalarge;
  const float lowerres = smallest[index] + frac*deltasmall;
  
  min_sigma_width_kev = min_width_multiple * lowerres;
  max_sigma_width_kev = max_width_multple * largerres;
  
  // For really nice HPGe or micro-calorimeters, the resolution may be even better than
  //  expected from the above, so we'll also check the spectrum an allow the FWHM
  //  go down to 2.12 channels/sigma (5 channels/FWHM).
  //  At lower energies (~100 keV), nice HPGe are easily ~1.5 channels/sigma
  //  so if anything we could probably cut this value down to 1.5 safely...
  const double min_nchannel_sigma = 2.12;  //2.12*2.355=4.99
  if( highres && meas && (meas->num_gamma_channels() > (4096+2)) )
  {
    shared_ptr<const SpecUtils::EnergyCalibration> cal = meas->energy_calibration();
    if( cal && cal->valid() )
    {
      try
      {
        const double highchannel = cal->channel_for_energy( energy + 0.5*min_sigma_width_kev );
        const double lowchannel = cal->channel_for_energy( energy - 0.5*min_sigma_width_kev );
        const double nchannel_sigma = highchannel - lowchannel;
        if( nchannel_sigma > min_nchannel_sigma )
          min_sigma_width_kev = min_sigma_width_kev * min_nchannel_sigma / nchannel_sigma;
      }catch( std::exception & )
      {
        // probably wont ever get here
      }//try / catch
    }//if( valid energy cal )
  }//if( highres )
}//void expected_peak_width_limits(...)



//combine_peaks_to_roi: throws exception on error
void combine_peaks_to_roi( PeakShrdVec &coFitPeaks,
                          double &roiLower,
                          double &roiUpper,
                          bool &lowstatregion,
                          const std::shared_ptr<const Measurement> &dataH,
                          const PeakShrdVec &inpeaks,
                          const double mean0,
                          const double sigma0,
                          const double area0,
                          const double pixelPerKev,
                          const bool isHPGe )
{
  typedef std::shared_ptr<const PeakDef> PeakDefShrdPtr;
  
  const double low_res_avrg_min_uncert_lowstat = 0.15;
  const double lowres_overlap_min_frac_to_combine = 0.2;
  const double lowres_nsigma_apart_to_combine = 5.0;
  const double lowres_max_nsigma_apart_to_combine = 10.0;
  
  assert( dataH );
  coFitPeaks.clear();
  roiLower = roiUpper = -1.0;
  
  const size_t nchannels = dataH->num_gamma_channels();
  double minEnergy = mean0 - 2.0*sigma0 - 20.0/pixelPerKev;
  double maxEnergy = mean0 + 2.0*sigma0 + 20.0/pixelPerKev;
  
  double roiLowerFeature, roiUpperFeature;
  
  {//begin codeblock to determine ROI limits
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "Intial mean=" << mean0 << ", sigma0=" << sigma0 << "\n";
#endif
    
    PeakDef dummypeak( mean0, sigma0, area0 );
    findROIEnergyLimits( roiLowerFeature, roiUpperFeature, dummypeak, dataH, isHPGe );

    roiLower = roiLowerFeature;
    roiUpper = roiUpperFeature;

    
    const size_t lchannel = dataH->find_gamma_channel( mean0-sigma0 );
    const size_t rchannel = dataH->find_gamma_channel( mean0+sigma0 );
    const double ndata = dataH->gamma_channels_sum( lchannel, rchannel );
    const double ndata_per_channel = ndata / (rchannel - lchannel + 1);
    const double avrg_uncert = 1.0 / sqrt(ndata_per_channel);
    
    lowstatregion = (avrg_uncert >= low_res_avrg_min_uncert_lowstat);
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "avrg_uncert=" << avrg_uncert << "\n"
                   << "Goes to roiLower=" << roiLower << ", roiUpper=" << roiUpper << "\n";
#endif

    if( !isHPGe )
    {
      try
      {
        find_roi_for_2nd_deriv_candidate( roiLower, roiUpper, mean0, dataH, isHPGe );
      }catch( std::exception &e )
      {
        //cerr << "find_roi_for_2nd_deriv_candidate failed: " << e.what() << endl;
        throw runtime_error( string("find_roi_for_2nd_deriv_candidate failed: ") + e.what() );
      }
    }//if( !isHPGe )
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "after Goes to roiLower=" << roiLower << ", roiUpper=" << roiUpper << ", sigma0=" << sigma0 << "\n";
#endif
  }//end codeblock to determine ROI limits
  
  
  PeakShrdVec nearPeaks = peaksTouchingRange( roiLower, roiUpper, inpeaks );
  
  //What we should actually do here is branch (to other functions) based on if
  //  we are going to be fitting for multiple peaks, or a single peak.  This way
  //  can can attempt to determine if we get a better result together or apart.
  
  for( const PeakDefShrdPtr &peak : nearPeaks )
  {
    if( peak->type() != PeakDef::GaussianDefined )
      continue;
    
    double overlapStart, overlapEnd;
    double peakmin = peak->lowerX();
    double peakmax = peak->upperX();
    
    //Doing this can cause the regions to no longer overlap
    //    if( !isHPGe )
    //      try{ find_roi_for_2nd_deriv_candidate( peakmin, peakmax, peak->mean(), dataH ); }catch(...){ }
    
    if( roiUpper>=peakmin && roiUpper<=peakmax )
    {
      overlapStart = peakmin;
      overlapEnd   = roiUpper;
    }else if( roiLower>=peakmin && roiUpper>=peakmax )
    {
      overlapStart = roiLower;
      overlapEnd   = peakmax;
    }else if( roiLower<=peakmin && roiUpper>=peakmax )
    {
      overlapStart = peakmin;
      overlapEnd   = peakmax;
    }else if( roiLower>=peakmin && roiUpper<=peakmax )
    {
      overlapStart = roiLower;
      overlapEnd   = roiUpper;
    }else
      throw runtime_error( "searchForPeakFromUser: logic error in range check" );
    
    const double nsigmaapart = fabs(peak->mean() - mean0) / peak->sigma();
    const double overlapAmount = overlapEnd - overlapStart;
    const double overlapNewFrac = overlapAmount / (roiUpper-roiLower);
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    const double overlapExistFrac = overlapAmount / (peakmax-peakmin);
    
    DebugLog(cerr) << "nsigmaapart=" << nsigmaapart << " for " << peak->mean() << " vs " << mean0 << "\n"
    << "overlapStart=" << overlapStart << ", overlapEnd=" << overlapEnd
    << ", overlapAmount=" << overlapAmount
    << ", overlapExistFrac=" << overlapExistFrac
    << ", overlapNewFrac=" << overlapNewFrac
    << ", npeaksigma=" << (fabs(peak->mean() - mean0)/peak->sigma())
    << "\n";
#endif
    
    if( isHPGe )
    {
      //numbers below arent based on much, as ov yet
      const double highres_nsigma_apart_combine = 8.5;
      const double highres_nsigma_dont_combine = 10.0;
      const double highres_overlap_frac_to_combine = 0.35;
      
      const double dist = fabs(peak->mean() - mean0);
      const bool isnear = (((1.0+dist)/peak->sigma()) < highres_nsigma_apart_combine);
      const bool isfar = ((dist/peak->sigma()) > highres_nsigma_dont_combine);
      const bool isoverlapping = (overlapNewFrac > highres_overlap_frac_to_combine);
      
      
      if( isnear || (isoverlapping && !isfar) )
      {
        coFitPeaks.push_back( peak );
      
        roiUpper  = std::max( roiUpper, peak->upperX() );
        roiLower  = std::min( roiLower, peak->lowerX() );
        if( mean0 < peak->mean() )
          maxEnergy = std::min( maxEnergy, peak->mean() );
        else
          minEnergy = std::max( minEnergy, peak->mean() );
      }
    }else
    {
      if( (overlapNewFrac > lowres_overlap_min_frac_to_combine
           || nsigmaapart < lowres_nsigma_apart_to_combine )
         && (nsigmaapart < lowres_max_nsigma_apart_to_combine) )
      {
        coFitPeaks.push_back( peak );
        
        roiUpper  = std::max( roiUpper, peak->upperX() );
        roiLower  = std::min( roiLower, peak->lowerX() );
        
        if( mean0 < peak->mean() )
          maxEnergy = std::min( maxEnergy, peak->mean() );
        else
          minEnergy = std::max( minEnergy, peak->mean() );
      }
    }//if( isHPGe ) / else
  }//for( nearPeaks )
  
  //Should handle special case where the ROI wont be combined by the above, but
  //  there is actually a ROI with peaks on both sides of the current one.
  
  
  //Now we have to grab any peaks that share a ROI with any of the peaks we added
  for( const PeakDefShrdPtr &in : inpeaks )
  {
    if( std::count(coFitPeaks.begin(), coFitPeaks.end(), in) )
      continue;
    
    if( (!in->continuum()->isPolynomial())
       || (in->type() != PeakDef::GaussianDefined) )
      continue;
    
    for( PeakDefShrdPtr peak : coFitPeaks )
    {
      if( peak->continuum() == in->continuum() )
      {
        coFitPeaks.push_back( in );
        
        roiLower = std::min( roiLower, in->lowerX() );
        roiUpper = std::max( roiUpper, in->upperX() );
        
        if( mean0 > in->mean() )
        {
          roiLower  = std::min( roiLower, in->lowerX() );
          minEnergy = std::max( minEnergy, in->mean() );
        }else
        {
          roiUpper  = std::max( roiUpper, in->upperX() );
          maxEnergy = std::min( maxEnergy, in->mean() );
        }
        break;
      }//if( already->continuum() == in->continuum() )
    }//for( PeakDefShrdPtr already : coFitPeaks )
  }//for( const PeakDefShrdPtr &peak, inpeaks )
  
  std::sort( coFitPeaks.begin(), coFitPeaks.end(),
            &PeakDef::lessThanByMeanShrdPtr );
  
  //Lets check to see if we should re-estimate the ROI
  if( coFitPeaks.size() )
  {
    const bool newIsLowest = (mean0 < coFitPeaks.front()->mean());
    const bool newIsHighest = (mean0 > coFitPeaks.back()->mean());
    
    if( newIsLowest || newIsHighest )
    {
      double thisROILower, thisROIUpper;
      PeakDef dummy( mean0, sigma0, 100.0 );
      findROIEnergyLimits( thisROILower, thisROIUpper, dummy, dataH, isHPGe );
      
      roiLower = std::min( roiLower, thisROILower );
      roiUpper = std::max( roiUpper, thisROIUpper );

      if( newIsLowest )
        dummy = PeakDef( coFitPeaks.back()->mean(), std::max(coFitPeaks.back()->sigma(),sigma0), coFitPeaks.back()->amplitude() );
      else
        dummy = PeakDef( coFitPeaks.front()->mean(), std::max(coFitPeaks.front()->sigma(),sigma0), coFitPeaks.front()->amplitude() );
      
      findROIEnergyLimits( thisROILower, thisROIUpper, dummy, dataH, isHPGe );
      
      roiLower = std::min( roiLower, thisROILower );
      roiUpper = std::max( roiUpper, thisROIUpper );

      if( !isHPGe )
      {
        try
        {
          double testLeftLower, testLeftUpper, testRightLower, testRightUpper;
          find_roi_for_2nd_deriv_candidate( testLeftLower, testLeftUpper,
                                           coFitPeaks.front()->mean(), dataH, isHPGe );
          find_roi_for_2nd_deriv_candidate( testRightLower, testRightUpper,
                                           coFitPeaks.back()->mean(), dataH, isHPGe );
          roiLower = std::min( roiLower, testLeftLower );
          roiUpper = std::max( roiUpper, testRightUpper );
        }catch( std::exception & )
        {
          cerr << "Failed to find candidate peak for multipeak region - continuuing on" << endl;
        }
      }//if( !isHPGe )
    }//if( newIsLowest || newIsHighest )
  }//if( we should re-estimate the ROI )
}//void combine_peaks_to_roi(...)


void get_candidate_peak_estimates_for_user_click(
                                 double &sigma0, double &mean0, double &area0,
                                 const double x,
                                 const double pixelPerKev,
                                 const std::shared_ptr<const Measurement> &dataH,
                                 const bool isHPGe,
                                 const PeakShrdVec &inpeaks )
{
  typedef std::shared_ptr<PeakDef> PeakPtr;
  typedef std::shared_ptr<const PeakDef> PeakConstPtr;
  
  const double lower_energy_mult = 0.2;
  const double upper_energy_mult = 0.2;
  
  const size_t nchannels = dataH->num_gamma_channels();
  
  const size_t midbin = dataH->find_gamma_channel( x );
  
  const double lower_chan_sub = lower_energy_mult*nchannels;
  assert( lower_chan_sub > 0 );
  size_t lowchannel = ((lower_chan_sub < midbin) ? (midbin - static_cast<size_t>(std::round(lower_chan_sub))) : 0);
  
  const double upper_chan_sub = upper_energy_mult*nchannels;
  size_t highchannel = ((midbin + upper_chan_sub) >= nchannels) ? nchannels-1 : static_cast<size_t>(midbin + upper_chan_sub);
  
  float min_sigma_width_kev, max_sigma_width_kev;
  expected_peak_width_limits( x, isHPGe, dataH, min_sigma_width_kev, max_sigma_width_kev );

  
  
  const size_t lower_reasonable_channel = dataH->find_gamma_channel( x - 20*max_sigma_width_kev );
  const size_t upper_reasonable_channel = dataH->find_gamma_channel( x + 20*max_sigma_width_kev );
  
  //cout << "For x=" << x << " keV, {lowchannel=" << lowchannel
  //     << ", lower_reasonable_channel=" << lower_reasonable_channel
  //     << "}, {" << highchannel << ", upper_reasonable_channel=" << upper_reasonable_channel << "}"
  //     << endl;
  
  lowchannel = std::max(lowchannel, lower_reasonable_channel);
  highchannel = std::min( highchannel, upper_reasonable_channel );
  
  
  const vector<PeakPtr> candidates
       = secondDerivativePeakCanidatesWithROI( dataH, isHPGe, lowchannel, highchannel );
  

  float min_sigma, max_sigma;
  expected_peak_width_limits( x, isHPGe, dataH, min_sigma, max_sigma );

  sigma0 = 0.5*(min_sigma + max_sigma) * (isHPGe ? 0.20 : 0.25);  //expected_peak_width_limits multiplies max width by  4 for isHPGe, and 3 for lowres
  mean0 = x;
  area0 = 100.0;
  
  bool updatedSigmaFromPrev = false;
  
  if( !updatedSigmaFromPrev )
  {
    shared_ptr<const PeakDef> leftpeak, rightpeak;
    
    for( const shared_ptr<const PeakDef> &p : inpeaks )
    {
      if( !p->gausPeak() )
        continue;
      
      if( (p->mean() <= x) && (!leftpeak || fabs(p->mean() - x) < fabs(leftpeak->mean() - x)) )
        leftpeak = p;
      
      if( (p->mean() >= x) && (!rightpeak || fabs(p->mean() - x) < fabs(rightpeak->mean() - x)) )
        rightpeak = p;
    }
    
    if( leftpeak && rightpeak )
    {
      const double dist_between = rightpeak->mean() - leftpeak->mean();
      const double frac_between = ((x - leftpeak->mean()) / dist_between);
      const double sigma_diff = rightpeak->sigma() - leftpeak->sigma();
      
      updatedSigmaFromPrev = true;
      sigma0 = leftpeak->sigma() + frac_between*sigma_diff;
    }
  }//if( !updatedSigmaFromPrev )
  
  
  if( !updatedSigmaFromPrev )
  {
    for( const shared_ptr<const PeakDef> &p : inpeaks )
    {
      if( p->gausPeak() && (fabs(p->mean() - x) < (10*p->sigma())) ) //10 is arbitrary.
      {
        updatedSigmaFromPrev = true;
        sigma0 = p->sigma();
      }
    }//for( const std::shared_ptr<const PeakDef> &p : inpeaks )
  }//if( !updatedSigmaFromPrev )
  
  if( sigma0 <= 0.0 )
    sigma0 = 1.0;  //JIC, shouldnt ever happen
  
  if( candidates.size() )
  {
    map<double,PeakPtr> candidatesMap;
    for( const PeakPtr &p : candidates )
    {
      //make sure candidate isnt within 0.75 sigma of existing peaks
      bool nearmean = false;
      for( const PeakConstPtr &in : inpeaks )
      {
        if( in->gausPeak() && fabs((in->mean()-p->mean())/std::min(in->sigma(),p->sigma())) < 0.75 )
        {
          nearmean = true;
          break;
        }
      }
      
      if( !nearmean )
      {
        for( auto dataDefPeak : inpeaks )
        {
          if( !dataDefPeak || dataDefPeak->gausPeak() )
            continue;
          
          //Candidate mean shouldnt be within a data defined peak ROI, and
          //  the candidate shouldnt span over the data defined peak.
          if( ((p->mean() > dataDefPeak->lowerX()) && (p->mean() < dataDefPeak->upperX()))
              || ((p->lowerX() < dataDefPeak->lowerX()) && (p->upperX() > dataDefPeak->upperX())) )
          {
            nearmean = true;
            break;
          }
      
          if( (p->lowerX() > dataDefPeak->lowerX()) && (p->lowerX() < dataDefPeak->upperX()) )
            p->continuum()->setRange( dataDefPeak->upperX(), p->upperX() );
          
          if( (p->upperX() > dataDefPeak->lowerX()) && (p->upperX() < dataDefPeak->upperX()) )
            p->continuum()->setRange( p->lowerX(), dataDefPeak->lowerX() );
        }//for( auto dataDefPeak : inpeaks )
      }//if( !nearmean )
      
      if( !nearmean )
        candidatesMap[fabs(p->mean()-x)] = p;
    }//for( const PeakDef &p : candidates )
    
    if( candidatesMap.size() )
    {
      const PeakDef &peak = *(candidatesMap.begin()->second);
      const double pixelUncert = 1.5*peak.sigma()*pixelPerKev + 20.0;
      const double pixelDelta = fabs(peak.mean()-x)*pixelPerKev;
      
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "For peak at " << peak.mean() << ", pixelDelta=" << pixelDelta
      << ", pixelUncert=" << pixelUncert << ", pixelPerKev=" << pixelPerKev
      << ", peak.sigma=" << (peak.gausPeak() ? peak.sigma() : 0.25*peak.roiWidth()) << "\n";
#endif
      if( (pixelDelta < pixelUncert) && pixelDelta < 75.0 )  //The 75 is arbitrary
      {
        sigma0 = peak.sigma();
        mean0 = peak.mean();
        area0 = peak.amplitude();
      }//if( pixelDelta < pixelUncert )
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      else
        DebugLog(cerr) << "Using dumb initial starting value guesses\n";
#endif
    }//if( candidatesMap.size() )
  }//if( candidates.size() )
  
}//get_candidate_peak_estimates_for_user_click(...)


//fit_peak_for_user_click(): upon error, 'results' will be empty.
//  Gives fit peaks in 'results', and the chi2/DOF of the fit in 'chiDof'.
//  All fit peaks will share the same continuum.
//  The initial guess of mean/sigma for the new peak being fit for must
//  be specified (area is currenly ignored), and all other peaks that will share
//  the ROI/continuum are specified in 'coFitPeaks'.
//  'lowerEnergies' and 'upperEnergies' specify the upper and lower bounds of
//  the ROI to try to fit for; in case more than one set of ROI limits is
//  supplied, the ROI range with the best chi2 will be used;  these input
//  vectors must be non-empty, and the same size.
void fit_peak_for_user_click( PeakShrdVec &results,
                              double &chi2Dof,
                              const std::shared_ptr<const Measurement> &dataH,
                              PeakShrdVec coFitPeaks,
                              const double mean0, const double sigma0,
                              const double area0,
                              const vector<double> &lowerEnergies,
                              const vector<double> &upperEnergies,
                             const bool isHPGe )
{
  typedef std::shared_ptr<const PeakDef> PeakDefShrdPtr;
  
  assert( !lowerEnergies.empty() );
  assert( lowerEnergies.size() == upperEnergies.size() );
  
  chi2Dof = DBL_MAX;
  results.clear();
  
  const size_t nchannels = dataH->num_gamma_channels();
  const size_t midbin = dataH->find_gamma_channel( mean0 );
  const float binwidth = dataH->gamma_channel_width( midbin );
  const size_t nFitPeaks = coFitPeaks.size() + 1;
  
  //The below should probably go off the number of bins in the ROI
  PeakContinuum::OffsetType offsetType;
  if( isHPGe )
    offsetType = (nFitPeaks < 3) ? PeakContinuum::Linear : PeakContinuum::Quadratic;
  else
    offsetType = (nFitPeaks < 2) ? PeakContinuum::Linear : PeakContinuum::Quadratic;
  
  for( size_t i = 0; i < coFitPeaks.size(); ++i )
    offsetType = std::max( offsetType, coFitPeaks[i]->continuum()->type() );
  
  std::shared_ptr<PeakDef> candidatepeak
                           = std::make_shared<PeakDef>(mean0, sigma0, area0);
  
  coFitPeaks.push_back( candidatepeak );
  std::sort( coFitPeaks.begin(), coFitPeaks.end(),
            &PeakDef::lessThanByMeanShrdPtr );
  
  for( size_t i = 0; i < lowerEnergies.size(); ++i )
  {
    const float lenergy = lowerEnergies[i];
    const float uenergy = upperEnergies[i];
    
    try
    {
      if( lenergy > uenergy )
        continue;
      
      const double range = (uenergy - lenergy);
      
      double minsigma = binwidth;
      double maxsigma = 0.5*range;
      
      ROOT::Minuit2::MnUserParameters params;
      
      for( size_t i = 0; i < coFitPeaks.size(); ++i )
      {
        char name[64];
        snprintf( name, sizeof(name), "Mean%i", static_cast<int>(i) );
        const double mean = coFitPeaks[i]->mean();
        const double sigma = coFitPeaks[i]->sigma();
        if( !coFitPeaks[i]->fitFor(PeakDef::Mean) )
        {
          params.Add( name, mean );
        }else
        {
          if( coFitPeaks[i] == candidatepeak )
            params.Add( name, mean, 0.1*sigma, mean-sigma, mean+sigma );
          else
            params.Add( name, mean, 0.1*range, mean-0.25*sigma, mean+0.25*sigma );
        }
        
        
        if( !isHPGe )
        {
//          cout << "Testing setting peak resolution limits based on expected_lowres_peak_width_limits" << endl;
          float lowersigma, uppersigma;
          expected_peak_width_limits( mean, false, dataH, lowersigma, uppersigma );
          if( !i )
            minsigma = lowersigma;
          if( i == (coFitPeaks.size()-1) )
            maxsigma = 1.33*uppersigma;
        }
      }//for( const PeakDefShrdPtr &peak : coFitPeaks )
      
      params.Add( "Sigma0", sigma0, 0.1*sigma0, minsigma, maxsigma );
      if( coFitPeaks.size() > 1 )
        params.Add( "SigmaFcn", 0.0, 0.001, -0.10, 0.10 );
            
      const PeakDef::SkewType skew_type = LinearProblemSubSolveChi2Fcn::skewTypeFromPrevPeaks( coFitPeaks );
      LinearProblemSubSolveChi2Fcn::addSkewParameters( params, skew_type, coFitPeaks );
      
      LinearProblemSubSolveChi2Fcn chi2Fcn( coFitPeaks, dataH, offsetType, skew_type,
                                           lenergy, uenergy );
      
      ROOT::Minuit2::MnUserParameterState inputParamState( params );
      ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
      
      ROOT::Minuit2::CombinedMinimizer fitter;
      unsigned int maxFcnCall = 0;
      const double tolerance = 0.01;
      
      ROOT::Minuit2::FunctionMinimum minimum
      = fitter.Minimize( chi2Fcn, params, strategy, maxFcnCall, tolerance );
      
      params = minimum.UserState().Parameters();
      
      const vector<double> pars = params.Params();
      const vector<double> errors = params.Errors();
      
      vector<PeakDef> fitpeaks;
      double firstFitChi2 = chi2Fcn.parametersToPeaks( fitpeaks, &pars[0], &errors[0] );
      double thisChi2Dof = firstFitChi2 / chi2Fcn.dof();
      
      for( const PeakDef peak : fitpeaks )
      {
        const double mean = peak.mean();
        const double sigma = peak.sigma();
        const double data = dataH->gamma_integral( mean-0.5*sigma, mean+0.5*sigma );
        const double area = peak.gauss_integral( mean-0.5*sigma, mean+0.5*sigma );
        
        if( area < 5.0 || area < sqrt(data) )
        {
          thisChi2Dof = DBL_MAX;
          break;
        }//if( area < 5.0 || area < sqrt(data) )
      }//for( const PeakDefShrdPtr peak : answer.first )
      
      
      if( thisChi2Dof < chi2Dof )
      {
        chi2Dof = thisChi2Dof;
        results.clear();
        for( const PeakDef &p : fitpeaks )
          results.push_back( std::make_shared<PeakDef>(p) );
      }//if( thisChi2Dof < chi2Dof )
    }catch( std::exception &e )
    {
      cout << "Caught: " << e.what() << endl;
    }//try / catch
  }//for( const float lenergy : lowerEnergies )
}//void fit_peak_for_user_click(...)


//lowres_shrink_roi(): all peaks in 'inpeaks' are assumed to share a ROI, and
//  have the correct definition of chi2/DOF.
//  This function has not been tested well, and could use further development.
//If there is an error 'inpeaks' will be returned.
PeakShrdVec lowres_shrink_roi( const PeakShrdVec &inpeaks,
                               const std::shared_ptr<const Measurement> &dataH,
                               const bool lowstatregion,
                               const bool automated )
{
  PeakShrdVec answer;
  
//  const double lowres_roi_shrink_min_nsigma_lowstat = automated ? 3.5 : 2.0;
//  const double lowres_roi_shrink_min_nsigma_highstat = 2.5;
  const double lowres_roi_rel_len_ratio = 1.25;
  
  if( inpeaks.empty() )
    return answer;
  
  try
  {
    vector<PeakDef> fitpeaks;
    for( size_t i = 0; i < inpeaks.size(); ++i )
      fitpeaks.push_back( *inpeaks[i] );
    
    double roiLower = inpeaks[0]->continuum()->lowerEnergy();
    double roiUpper = inpeaks[0]->continuum()->upperEnergy();
    
    //  now go through and try to ahrink ROI by taking off bins until the end bin is close
    //  to the average chi2/bin of the fit, or when datahits zero or when gaussian
    //  contribution becomes significant
    const size_t initialFirstChannel = dataH->find_gamma_channel( roiLower );
    const size_t initialLastChannel = dataH->find_gamma_channel( roiUpper );
      
    size_t finalFirstChannel = initialFirstChannel;
    size_t finalLastChannel = initialLastChannel;
      
    double finalUpperE = roiUpper;
    double finalLowerE = roiLower;
    const double origFitChi2Dof = inpeaks[0]->chi2dof();
    double chi2Dof = origFitChi2Dof;
    
    //Try to shrink the upper edge of ROI
    //We could now look at the the regions above ~4 sigma away, and see if
    //  their chi2/dof is just really good, in which case we should consider
    //  removing.  But for right now (20141203) the above appears to be
    //  working mostly pretty well, see
    //  W187_GR135P(NaI)_52sec_0.5cmW.chn_20140724T124116.n42
    //  (upper ROIO range for peak at 777 keV) for a case that doesnt work so
    //  well currently (not bad, just ROI is larger than necassarry).
    const double lower_sigma = fitpeaks.front().sigma();
    const double upper_sigma = fitpeaks.back().sigma();
    const double lower_start = (fitpeaks.front().mean() - 2.0*lower_sigma);
    const double upper_start = (fitpeaks.back().mean() + 2.0*upper_sigma);
    
    double lower_extent = std::max( 2.0*lower_sigma, lower_start - finalLowerE );
    double upper_extent = std::max( 2.0*upper_sigma, finalUpperE - upper_start );
      
    if( upper_extent > lowres_roi_rel_len_ratio*lower_extent )
    {
      const double startx = upper_start + lower_extent;
      const double endx = finalUpperE;
      const double tailChi2Dof = evaluate_chi2dof_for_range( fitpeaks, dataH, startx, endx );
      const double generalChi2Dof = evaluate_chi2dof_for_range( fitpeaks, dataH, finalLowerE, startx );
      
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cout) << "Above tail checked " << startx << " to " << endx << "\n"
           << "tailChi2Dof=" << tailChi2Dof
           << ", generalChi2Dof=" << generalChi2Dof
           << ", origFitChi2Dof=" << origFitChi2Dof << "\n";
#endif
      
      if( tailChi2Dof < 0.95*generalChi2Dof )
      {
        double newChi2Dof;
        vector<PeakDef> newfitpeaks;
        refit_for_new_roi( inpeaks, dataH, finalLowerE, startx,
                           newChi2Dof, newfitpeaks );
          
        const double totalOldlen  = finalUpperE - finalLowerE;
        const double totalNewLen  = startx - finalLowerE;
        const double totalRemovedFrac = (totalOldlen - totalNewLen) / totalOldlen;
        double oldChi2Body = (chi2Dof - totalRemovedFrac*tailChi2Dof) /( 1.0 - totalRemovedFrac );
        
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cout) << "newChi2Dof=" << newChi2Dof << ", oldchi2dof=" << chi2Dof
             << ", oldChi2Body=" << oldChi2Body << "\n";
#endif
        
        //1.25 is arbitrary
        if( newChi2Dof < 1.25*oldChi2Body )
        {
          chi2Dof = newChi2Dof;
          fitpeaks = newfitpeaks;
          finalLastChannel = dataH->find_gamma_channel( startx );
          finalUpperE = startx;
        }
      }//if( tailChi2Dof < 0.95*generalChi2Dof )
        
    }else if( lower_extent > lowres_roi_rel_len_ratio*upper_extent )
    {
      const double startx = finalLowerE;
      const double endx = finalLowerE + upper_extent;
      const double tailChi2Dof = evaluate_chi2dof_for_range( fitpeaks, dataH, startx, endx );
      const double generalChi2Dof = evaluate_chi2dof_for_range( fitpeaks, dataH, endx, finalUpperE );
      
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cout) << "Below tail checked " << startx << " to " << endx << "\n"
           << "tailChi2Dof=" << tailChi2Dof
           << ", generalChi2Dof=" << generalChi2Dof
           << ", chi2Dof=" << chi2Dof << "\n";
#endif
      
      //0.95 is arbitrary
      if( tailChi2Dof < 0.95*generalChi2Dof )
      {
        double newChi2Dof;
        vector<PeakDef> newfitpeaks;
        refit_for_new_roi( inpeaks, dataH, endx, finalUpperE,
                            newChi2Dof, newfitpeaks );
          
        const double totalOldlen  = finalUpperE - finalLowerE;
        const double totalNewLen  = finalUpperE - endx;
        const double totalRemovedFrac = (totalOldlen - totalNewLen) / totalOldlen;
        double oldChi2Body = (chi2Dof - totalRemovedFrac*tailChi2Dof) /( 1.0 - totalRemovedFrac );
          
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cout) << "newChi2Dof=" << newChi2Dof << ", oldchi2dof=" << chi2Dof << ", oldChi2Body=" << oldChi2Body << "\n";
#endif
        
        //1.25 is arbitrary
        if( newChi2Dof < 1.25*oldChi2Body )
        {
          chi2Dof = newChi2Dof;
          fitpeaks = newfitpeaks;
          finalFirstChannel = dataH->find_gamma_channel( endx );
          finalLowerE = endx;
        }
      }//if( tailChi2Dof < 0.95*generalChi2Dof )
    }//if( upper_extent > 1.5*lower_extent ) / else
    
    if( fitpeaks.size() )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cout) << "Stopping shrinking upper end at channel " << finalLastChannel << " {"
           << dataH->gamma_channel_lower( finalLastChannel )
           << " kev}, (started at " << initialLastChannel << " {"
        << dataH->gamma_channel_lower(initialLastChannel) << " kev}), newChi2Dof="
        << chi2Dof << ", origFitChi2Dof=" << origFitChi2Dof
        << ", preAmp=" << inpeaks.back()->amplitude()
        << ", not Amp=" << fitpeaks.back().amplitude() << " keV"
        << "\n";
#endif
      
      if( (finalLastChannel != initialLastChannel)
            || (finalFirstChannel != initialFirstChannel) )
      {
        double newChi2Dof;
        vector<PeakDef> newfitpeaks;
          
        PeakShrdVec peaksToReFit;
        for( const PeakDef &p : fitpeaks )
          peaksToReFit.push_back( std::make_shared<PeakDef>(p) );
          
        const vector<double> originalFitPars, originalFitErrors;
          
        refit_for_new_roi( peaksToReFit, dataH, finalLowerE, finalUpperE,
                            newChi2Dof, newfitpeaks );
          
          
        if( (newChi2Dof < chi2Dof) && (newfitpeaks.size() > 0) )
        {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          cerr << "Refitign peak dropped chi2/dof from " << chi2Dof << " to " << newChi2Dof << endl;
#endif
          chi2Dof = newChi2Dof;
          fitpeaks = newfitpeaks;
        }else
        {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          cerr << "Refit of peaks didnt actually help! newChi2Dof=" << newChi2Dof << " vs chi2Dof=" << chi2Dof << endl;
#endif
        }
        
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cout) << "Final PeakRoiLower=" << fitpeaks[0].continuum()->lowerEnergy()
             << ", PeakRoiUpper=" << fitpeaks[0].continuum()->upperEnergy() << "\n";
#endif
      }//if( finalLastBin != upBin )
        
      answer.clear();
      for( size_t i = 0; i < fitpeaks.size(); ++i )
        answer.push_back( std::make_shared<PeakDef>(fitpeaks[i]) );
    }//if( fitpeaks.size() )
  }catch( std::exception & )
  {
    answer.clear();
  }
  
  if( answer.empty() )
    return inpeaks;
  
  return answer;
}//PeakShrdVec lowres_shrink_roi( PeakShrdVec initialfitpeaks )



PeakShrdVec highres_shrink_roi( const PeakShrdVec &inpeaks,
                                const std::shared_ptr<const Measurement> &dataH,
                                const bool lowstatregion,
                                const bool automated )
{
  PeakShrdVec answer;
  
  const double nsigma_start_shrinking = 10.0;
  const double nsigma_shrink_to = 8.0;
  
  if( inpeaks.empty() )
    return answer;
  
  try
  {
    vector<PeakDef> fitpeaks;
    for( size_t i = 0; i < inpeaks.size(); ++i )
      fitpeaks.push_back( *inpeaks[i] );
    
    const double origRoiLower = inpeaks[0]->continuum()->lowerEnergy();
    const double origRoiUpper = inpeaks[0]->continuum()->upperEnergy();
    
    const size_t initialFirstChannel = dataH->find_gamma_channel( origRoiLower );
    const size_t initialLastChannel = dataH->find_gamma_channel( origRoiUpper );
    
    size_t finalFirstChannel = initialFirstChannel;
    size_t finalLastChannel = initialLastChannel;
    
    double finalUpperE = origRoiUpper;
    double finalLowerE = origRoiLower;
    const double origFitChi2Dof = inpeaks[0]->chi2dof();
    double chi2Dof = origFitChi2Dof;
    
    //Try to shrink the upper edge of ROI
    //We could now look at the the regions above ~4 sigma away, and see if
    //  their chi2/dof is just really good, in which case we should consider
    //  removing.  But for right now (20141203) the above appears to be
    //  working mostly pretty well, see
    //  W187_GR135P(NaI)_52sec_0.5cmW.chn_20140724T124116.n42
    //  (upper ROIO range for peak at 777 keV) for a case that doesnt work so
    //  well currently (not bad, just ROI is larger than necassarry).
    const double lower_mean = fitpeaks.front().mean();
    const double upper_mean = fitpeaks.back().mean();
    const double lower_sigma = fitpeaks.front().sigma();
    const double upper_sigma = fitpeaks.back().sigma();
    
//    const double tailChi2Dof = evaluate_chi2dof_for_range( fitpeaks, dataH, startx, endx );
//    const double generalChi2Dof = evaluate_chi2dof_for_range( fitpeaks, dataH, finalLowerE, startx );
    
    
    const double nLowOriginal = (lower_mean - origRoiLower) / lower_sigma;
    const double nHighOriginal = (origRoiUpper - upper_mean) / upper_sigma;
    
    const bool tryShrinkLow = (nLowOriginal > nsigma_start_shrinking);
    const bool tryShrinkHigh = (nHighOriginal > nsigma_start_shrinking);
    
    if( tryShrinkLow || tryShrinkHigh )
    {
      //const double origChi2dof = inpeaks[0]->chi2dof();
      
//      const double extraTailChi2Dof = evaluate_chi2dof_for_range( fitpeaks, dataH, origRoiLower, lower_mean-3.0*nsigmatry );
      const double test_lower_roi = tryShrinkLow ? lower_mean-nsigma_shrink_to*lower_sigma : origRoiLower;
      const double test_upper_roi = tryShrinkHigh ? upper_mean+nsigma_shrink_to*upper_sigma : origRoiUpper;
      
      double testChi2Dof;
      vector<PeakDef> newfitpeaks;
      refit_for_new_roi( inpeaks, dataH, test_lower_roi, test_upper_roi, testChi2Dof, newfitpeaks );
      const double origChi2ForNewRange = evaluate_chi2dof_for_range( fitpeaks, dataH, test_lower_roi, test_upper_roi );
      
//      if( fabs(lower_mean-144.176) < 1.0 )
//        cout << endl;
      
      if( newfitpeaks.size() && (testChi2Dof < (1.1*origChi2ForNewRange+0.2)) )
      {
        const double delta_lower_mean = fabs(lower_mean - newfitpeaks.front().mean())/lower_sigma;
        const double delta_upper_mean = fabs(upper_mean - newfitpeaks.back().mean())/upper_sigma;
        
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cout) << "Shrinking from {" <<origRoiLower
        << ", " << origRoiUpper << "} to {" << test_lower_roi << ", "
        << test_upper_roi << "} made chi2dof go from " << origFitChi2Dof << " to "
        << testChi2Dof << ", delta_lower_mean=" << delta_lower_mean
        << ", delta_upper_mean=" << delta_upper_mean << "\n";
#endif

        if( delta_lower_mean < 0.05 && delta_upper_mean < 0.05 )
        {
          try
          {
            bool lowisgood = true, highisgood = true;
            
            const vector<float> &energies = *dataH->gamma_channel_energies();
            const vector<float> &contents = *dataH->gamma_channel_contents();
            
            if( tryShrinkLow )
            {
              const size_t start_channel = initialFirstChannel;
              const size_t end_channel = dataH->find_gamma_channel( lower_mean - 5.0*lower_sigma );
              
              std::vector<double> poly_coeffs, coeff_uncerts;
              const size_t nchannelsrange = (1 + end_channel - start_channel);
              const float *x = &(energies[0]) + start_channel;
              const float *y = &(contents[0]) + start_channel;
              const double linechi2 = fit_to_polynomial( x, y, nchannelsrange,
                                                1, poly_coeffs, coeff_uncerts );
              const double linchi2dof = linechi2 / (nchannelsrange - 2);
              
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
              cout << lower_mean << ", Lowerside linchi2dof=" << linchi2dof << ", nLowOriginal=" << nLowOriginal << endl;
#endif
              
              lowisgood = (linchi2dof < 1.5);
            }//if( tryShrinkLow )
            
            if( tryShrinkHigh )
            {
              const size_t start_channel = dataH->find_gamma_channel( upper_mean + 5.0*upper_sigma );
              const size_t end_channel = initialLastChannel;
              
              std::vector<double> poly_coeffs, coeff_uncerts;
              const size_t nchannelsrange = (1 + end_channel - start_channel);
              const float *x = &(energies[0]) + start_channel;
              const float *y = &(contents[0]) + start_channel;
              const double linechi2 = fit_to_polynomial( x, y, nchannelsrange,
                                                        1, poly_coeffs, coeff_uncerts );
              const double linchi2dof = linechi2 / (nchannelsrange - 2);
              
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
              cout << upper_mean << ", Upperside linchi2dof=" << linchi2dof <<", nHighOriginal=" << nHighOriginal << endl;
#endif
              
              highisgood = (linchi2dof < 1.5);
            }//if( tryShrinkHigh )
          
            //fit a line from
            if( lowisgood && highisgood )
            {
              finalFirstChannel = dataH->find_gamma_channel( test_lower_roi );
              finalLastChannel = dataH->find_gamma_channel( test_upper_roi );
              finalLowerE = test_lower_roi;
              finalUpperE = test_upper_roi;
              fitpeaks = newfitpeaks;
            } else
            {
              //Check if its a super high stat region - if so, try to shrink the
              //  ROI even more
            }
          }catch(...)
          {
          }
        }//if( delta_lower_mean < 0.05 && delta_upper_mean < 0.05 )
      }
      
    }//if( nLowOriginal > nsigma_start_shrinking )
    
  
    
    if( fitpeaks.size() )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cout) << "Stopping shrinking upper end at channel " << finalLastChannel << " {"
      << dataH->gamma_channel_lower( finalLastChannel )
      << " kev}, (started at " << initialLastChannel << " {"
      << dataH->gamma_channel_lower(initialLastChannel) << " kev}), newChi2Dof="
      << chi2Dof << ", origFitChi2Dof=" << origFitChi2Dof
      << ", preAmp=" << inpeaks.back()->amplitude()
      << ", not Amp=" << fitpeaks.back().amplitude() << " keV"
      << "\n";
#endif
      
      if( (finalLastChannel != initialLastChannel)
         || (finalFirstChannel != initialFirstChannel) )
      {
        double newChi2Dof;
        vector<PeakDef> newfitpeaks;
        
        PeakShrdVec peaksToReFit;
        for( const PeakDef &p : fitpeaks )
          peaksToReFit.push_back( std::make_shared<PeakDef>(p) );
        
        const vector<double> originalFitPars, originalFitErrors;
        
        refit_for_new_roi( peaksToReFit, dataH, finalLowerE, finalUpperE,
                          newChi2Dof, newfitpeaks );
        
        vector<PeakDef> origpeaks;
        for( size_t i = 0; i < inpeaks.size(); ++i )
          origpeaks.push_back( *inpeaks[i] );
//        const double origChi2ForNewRange = evaluate_chi2dof_for_range( origpeaks, dataH, finalLowerE, finalUpperE );
        
        
        if( (newChi2Dof < (chi2Dof+0.2)) && (newfitpeaks.size() > 0) )
        {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          cerr << "Refitign peak dropped chi2/dof from " << chi2Dof << " to " << newChi2Dof << endl;
#endif
          chi2Dof = newChi2Dof;
          fitpeaks = newfitpeaks;
        }else
        {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          cerr << "Refit of peaks didnt actually help! newChi2Dof=" << newChi2Dof << " vs chi2Dof=" << chi2Dof << endl;
#endif
        }
        
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cout) << "Final PeakRoiLower=" << fitpeaks[0].continuum()->lowerEnergy()
        << ", PeakRoiUpper=" << fitpeaks[0].continuum()->upperEnergy() << "\n";
#endif
      }//if( finalLastBin != upBin )
      
      answer.clear();
      for( size_t i = 0; i < fitpeaks.size(); ++i )
        answer.push_back( std::make_shared<PeakDef>(fitpeaks[i]) );
    }//if( fitpeaks.size() )
  }catch( std::exception &e )
  {
    answer.clear();
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "Caught exception; clearing peaks: " << e.what() << "\n";
#endif
  }
  
  if( answer.empty() )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "Answer is empty, returning inpeaks.\n";
#endif
    
    return inpeaks;
  }
  
  return answer;
}//PeakShrdVec highres_shrink_roi( PeakShrdVec initialfitpeaks )




bool check_lowres_single_peak_fit( const std::shared_ptr<const PeakDef> peak,
                                   const std::shared_ptr<const Measurement> &dataH,
                                   const bool lowstatregion,
                                   const bool automated )
{
  const double lowres_max_chi2dof = automated ? 25.0 : 50.0;
  const double lowres_min_core_chi2dof_peak_improvment = automated ? 0.5 : 0.25;
  const double lowres_min_withinsigma_chi2dof_peak_improvment = automated ? 1.5 : 0.5;
  const double lowres_min_energy_require_chi2dof_cut = 100.0;
  const double lowres_max_nsigma_to_require_chi2dof_cut = 50;
  const double lowres_min_nsigma_peak = automated ? 3.0 : 1.5;
  const double lowres_roi_min_nsigma = 3.5;
  const double lowres_min_for_narrow_roi_nsigma_peak = automated ? 10.0 : 5.0;
  const double lowres_min_core_chi2dof_over_line_improvment = automated ? 0.8 : 0.5;
  const double lowres_min_core_line_chi2dof_lower_ratio = automated ? 0.5 : 0.25;
  const double lowres_bad_continuum_fit_multiple = automated ? 0.333 : 0.0;
  const double lowres_min_necassary_chi2_improv_over_line = automated ? 0.35 : 0.1;
  const double lowres_badcont_line_chi2dof = 4.0;
  const double lowres_badcont_corecontinuum_chi2dof = 10.0;
  const double lowres_badcont_cont_to_line_ratio = 3.5;
  const bool lowres_enforce_peak_width_limits = true;
  const double max_avrg_uncert_to_require_chi2_improv_over_line = 0.45;
  
  const double fwhm = peak->fwhm();
  const double mean = peak->mean();
  const double sigma = peak->sigma();
  const double core_start = std::max( mean - fwhm, peak->lowerX() );
  const double core_end = std::min( mean + fwhm, peak->upperX() );
      
  vector<PeakDef> fitpeaks( 1, *peak ), fitpeaksnoamp( 1, *peak );
  fitpeaksnoamp[0].setAmplitude( 0.0 );
  const double core_chi2dof = evaluate_chi2dof_for_range( fitpeaks,
                                                             dataH, core_start, core_end );
  const double no_peak_chi2dof = evaluate_chi2dof_for_range( fitpeaksnoamp,
                                                                dataH, core_start, core_end );
  
  if( (no_peak_chi2dof - core_chi2dof) < lowres_min_core_chi2dof_peak_improvment )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "check_lowres_single_peak_fit: Failed to fit a peak for the "
         << "change in chi2 without the peak being to small (chi2dof="
         << no_peak_chi2dof << " w/ continuum only, " << core_chi2dof
         << " w/ continuum and peak)" << "\n";
#endif
    return false;
  }//if( (no_peak_chi2dof - core_chi2dof) < 0.5 )
      
  const double withpeakchi2dof = evaluate_chi2dof_for_range( fitpeaks, dataH, mean-sigma, mean+sigma );
  const double withoutpeakschi2dof = evaluate_chi2dof_for_range( fitpeaksnoamp, dataH, mean-sigma, mean+sigma );
        
  //Isnt this just a duplicate of the above???
  if( !lowstatregion && withoutpeakschi2dof < (withpeakchi2dof+lowres_min_withinsigma_chi2dof_peak_improvment) )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak for not improving"
         << " the chi2/dof enough withpeakchi2dof=" << withpeakchi2dof
         << ", withoutpeakschi2dof=" << withoutpeakschi2dof << "\n";
#endif
    
    return false;
  }
  
  const double gausarea = peak->gauss_integral( mean-2.0*sigma, mean+2.0*sigma );
  const double dataarea = dataH->gamma_integral( mean-2.0*sigma, mean+2.0*sigma );
  const double nsigma = gausarea / sqrt( std::max(dataarea,1.0) );
  const double chi2Dof = peak->chi2dof();

#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  DebugLog(cout) << "mea " << mean << ", nsigma=" << nsigma << "\n";
#endif
  
  if( chi2Dof > lowres_max_chi2dof
      && mean > lowres_min_energy_require_chi2dof_cut
      && nsigma < lowres_max_nsigma_to_require_chi2dof_cut )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak for the chi2DOF"
         << " being to bad (" << chi2Dof << ")" << "\n";
#endif
    
    return false;
  }//if( chi2Dof > 25.0 )
  
  if( nsigma < lowres_min_nsigma_peak )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak for the gros stat"
         << " significance to low (" << nsigma << ")" << "\n";
#endif
    
    return false;
  }
  
  //Check if ROI extent kinda makes sense for this peak
  const double roi_width = peak->upperX() - peak->lowerX();
      
  if( ((roi_width/sigma) < lowres_roi_min_nsigma)
     && (nsigma < lowres_min_for_narrow_roi_nsigma_peak) )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak for the ROI"
         << " being to narrow for this peak (roiwidth=" << roi_width
         << ", sigma=" << sigma << "), with nsigma only equal to "
         << nsigma << "\n";
#endif
    
    return false;
  }
  
//The following test erroneously fails CZT spectra actual peaks
//  const double lowres_roi_max_nsigma = 20.0;
//  cout << mean << ", (roi_width/sigma)=" << (roi_width/sigma) << endl;
//  if( (roi_width/sigma) > lowres_roi_max_nsigma )
//  {
//    cerr << "check_lowres_single_peak_fit: Failed to fit a peak for the sigma"
//    << " being too small for ROI extent (roiwidth=" << roi_width
//    << ", sigma=" << sigma << ")" << endl;
//    return false;
//  }
    
  
  //Make sure the peak width is within bounds of what would be expected
  if( lowres_enforce_peak_width_limits )
  {
    float min_sigma, max_sigma;
    expected_peak_width_limits( mean, false, dataH, min_sigma, max_sigma );
    
    if( sigma < min_sigma || sigma > max_sigma )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak for the width"
      << " being out of expected range "
      << "({mean: " << mean << ", sigma: " << sigma << ", min_sigma: "
      << min_sigma << ", max_sigma: " << max_sigma << "})" << "\n";
#endif
      
      return false;
    }
  }//if( lowres_enforce_peak_width_limits )
  
  
  //check that we couldnt just draw a line and have the same chi2 as with the
  //  peak
  try
  {
    const vector<float> &energies = *dataH->gamma_channel_energies();
    const vector<float> &contents = *dataH->gamma_channel_contents();
    const size_t lowerchannel = dataH->find_gamma_channel( peak->lowerX() );
    const size_t upperchannel = dataH->find_gamma_channel( peak->upperX() );
        
    std::vector<double> poly_coeffs, coeff_uncerts;
    const size_t nchannelsrange = (1 + upperchannel - lowerchannel);
    fit_to_polynomial( &(energies[0]) + lowerchannel,
                       &(contents[0]) + lowerchannel,
                       nchannelsrange, 1, poly_coeffs, coeff_uncerts );
    
    const size_t lowercore = dataH->find_gamma_channel( core_start );
    const size_t uppercore = dataH->find_gamma_channel( core_end );
        
    double linechi2 = 0.0, avrguncert = 0.0;
    for( size_t channel = lowercore; channel <= uppercore; ++channel )
    {
      const double x = dataH->gamma_channel_lower( channel );
      const double y = dataH->gamma_channel_content( channel );
      const double y_pred = poly_coeffs[0] + poly_coeffs[1]*x;
      const double uncert = (y > MIN_CHANNEL_UNCERT ? std::sqrt(y) : 1.0);
      const double thichi2 = std::pow( (y_pred - y) / uncert, 2.0 );
      
      linechi2 += thichi2;
      avrguncert += 1.0 / uncert;
      
//      cerr << "\t{" << x << ", " << y << ", " << y_pred << ", " << thichi2
//           << ", " << linechi2/(channel-lowercore+1) << "}" << endl;
    }//for( int bin = 0; bin < nbin; ++bin )
    
    avrguncert /= (uppercore - lowercore + 1);
    const double linechi2dof = linechi2 / (uppercore - lowercore + 1);

#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "For " << mean << " core_chi2dof=" << core_chi2dof
         << ", linechi2dof=" << linechi2dof << ", nchannelsrange=" << nchannelsrange << "\n"
     << "avrguncert=" << avrguncert << "\n";
#endif
    
    const double n1sigma = mean - sigma;
    const double p1sigma = mean + sigma;
    const double peakarea = peak->gauss_integral( n1sigma, p1sigma );
    const double dataarea = dataH->gamma_integral( n1sigma, p1sigma );
        
    size_t lbin = dataH->find_gamma_channel( n1sigma );
    size_t rbin = dataH->find_gamma_channel( p1sigma );
    const size_t nbin_pm1sigma = (rbin - lbin + 1);
        
    const double peak_above_uncert = peakarea / sqrt( dataarea );
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "peak_above_uncert=" << peak_above_uncert << " over nbins=" << nbin_pm1sigma
        << "--->" << (peak_above_uncert/nbin_pm1sigma) << "\n"
        << "For " << mean
         << "\tcore_chi2dof=" << core_chi2dof << "\n"
         << "\tlinechi2dof=" << linechi2dof << "\n"
         << "\tno_peak_chi2dof=" << no_peak_chi2dof << "\n"
         << "\t(linechi2dof-core_chi2dof)=" << (linechi2dof - core_chi2dof) << "\n"
        << "\t(linechi2dof-core_chi2dof)/linechi2dof=" << (linechi2dof - core_chi2dof)/linechi2dof << "\n"
        << "\t(core_chi2dof/linechi2dof)=" << (core_chi2dof/linechi2dof) << "\n"
        << "\t(linechi2dof/core_chi2dof)=" << (linechi2dof/core_chi2dof) << ", (no_peak_chi2dof/linechi2dof)=" << (no_peak_chi2dof/linechi2dof) << "\n";
#endif
    
    //if the continuum only chi2 is much different then the line-only chi2,
    if( (linechi2dof - core_chi2dof) < lowres_min_core_chi2dof_over_line_improvment
         && ((core_chi2dof/linechi2dof) > lowres_min_core_line_chi2dof_lower_ratio)
         && ((peak_above_uncert/nbin_pm1sigma) > 0.85
               || (linechi2dof/core_chi2dof) < lowres_bad_continuum_fit_multiple*(no_peak_chi2dof/linechi2dof)) )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
           << "change in chi2 over a straight line being to small (core_chi2dof="
           << core_chi2dof << ", linechi2dof=" << linechi2dof << ")"
           << "\n";
#endif
      return false;
    }//if( (no_peak_chi2dof - core_chi2dof) < 0.5 )
  
    if( (linechi2dof-core_chi2dof) < lowres_min_necassary_chi2_improv_over_line
       && avrguncert < max_avrg_uncert_to_require_chi2_improv_over_line )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
           << "change in chi2 over a straight line being to small no matter what"
           << " (core_chi2dof=" << core_chi2dof << ", linechi2dof="
           << linechi2dof << ")" << "\n";
#endif
      return false;
    }

    //This next test is intended to elimnate peaks where the continuum is
    //  much below where it visually should be, but because of the gausians
    //  amplitude the overall chi2 is decent in the fit with the peak, but
    //  really this is due to having the incorrect continuum.
    const double lowres_badcont_upper_chi2diff = 0.75;
    if( linechi2dof < lowres_badcont_line_chi2dof
           && no_peak_chi2dof > lowres_badcont_corecontinuum_chi2dof
           && (no_peak_chi2dof/linechi2dof) > lowres_badcont_cont_to_line_ratio
           && (linechi2dof-core_chi2dof) < lowres_badcont_upper_chi2diff )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
           << "a line fit the region well, but the peaks continuum didnt"
           << " (core_chi2dof=" << core_chi2dof  << ", no_peak_chi2dof="
           << no_peak_chi2dof << ", linechi2dof=" << linechi2dof
           << ")" << "\n";
#endif
      
      return false;
    }
    
    const double lowres_line_is_close = 1.1;
    const double lowres_line_is_close_necassary_improvment = 0.1;
    if( linechi2dof < lowres_line_is_close
        && (linechi2dof-core_chi2dof) < lowres_line_is_close_necassary_improvment )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
      << "a line fit the region well, and adding a peak didnt improve fit enough"
      << " (linechi2dof=" << linechi2dof  << ", core_chi2dof="
      << core_chi2dof << ")" << "\n";
#endif
      
      return false;
    }

    if( avrguncert < 0.2 && linechi2dof < 1.2 )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
      << "a line fit the region well, and average uncertainty wast too high"
      << " (linechi2dof=" << linechi2dof  << ", avrguncert="
      << avrguncert << ")" << "\n";
#endif
      
      return false;
    }
    
    
    //make sure the continuum isnt changing faster than the peak
    // (being lazy and just integrating, rather than evaluating)
    const double ptipval = peak->gauss_integral( mean-0.1*sigma, mean+0.1*sigma );
    const double p2val = peak->gauss_integral( mean+1.9*sigma, mean+2.1*sigma );
    const double conttipval = peak->offset_integral( mean-0.1*sigma, mean+0.1*sigma, dataH );
    const double cont2val = peak->offset_integral( mean+1.9*sigma, mean+2.1*sigma, dataH );
    const double contdiff = conttipval - cont2val;
    const double peakdiff = ptipval - p2val;
    
    const double max_relative_continuum_slope = 1.25;
    const double minuncert_apply_cont_slope_test = 0.05;
    
    const double roi_lower = peak->lowerX();
    const double roi_upper = peak->upperX();
    
    const double below_roi_cont_area = peak->continuum()->offset_integral( roi_lower - sigma, roi_lower, dataH );
    const double above_roi_cont_area = peak->continuum()->offset_integral( roi_upper, roi_upper + sigma, dataH );
    const double below_roi_data_area = dataH->gamma_integral( roi_lower - sigma, roi_lower );
    const double above_roi_data_area = dataH->gamma_integral( roi_upper, roi_upper + sigma );
    
    const double below_extra_nsigma = (below_roi_data_area - below_roi_cont_area) / sqrt(std::max(1.0,below_roi_data_area));
    const double above_extra_nsigma = (above_roi_data_area - above_roi_cont_area) / sqrt(std::max(1.0,above_roi_data_area));
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "above_extra_nsigma=" << above_extra_nsigma
         << ", below_extra_nsigma=" << below_extra_nsigma
         << ", contdiff/peakdiff=" << (contdiff/peakdiff) << "\n";
#endif
    
    //We will only require a not-so-fast drop in data for low statistics areas
    //  (since high stat continuums can fall faster than the peak drops off),
    //  and areas where the continuum goes above the data on both sides of the
    //  peak (this later condition catches compton bumbs)
    if( fabs(contdiff/peakdiff) > max_relative_continuum_slope
        && (avrguncert > minuncert_apply_cont_slope_test
            || (above_extra_nsigma<-3.0 && below_extra_nsigma<-3.0)) )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the"
      << " continuum is steeper than the peak"
      << "({peakdrop: " << peakdiff << ", continuumdrop: " << contdiff
      << "})"
      << " ratio=" << (contdiff/peakdiff) << ", nsigma=" << nsigma
      << ", avrguncert=" << avrguncert << "\n";
#endif
      
      return false;
    }//If slope is too steep )
    
  }catch( std::exception &e )
  {
    cerr << "check_lowres_single_peak_fit(): caught " << e.what() << endl;
  }
  
  return true;
}//bool check_lowres_single_peak_fit(...)


enum PeakRejectionStatus
{
  AcceptPeak,
  RejectedCanTryAgain,
  RejectedDontTryAgain
};//enum PeakRejectionStatus

PeakRejectionStatus check_lowres_multi_peak_fit( const vector<std::shared_ptr<const PeakDef> > fitpeaks,
                                  const vector<std::shared_ptr<const PeakDef> > originalpeaks,
                                  const std::shared_ptr<const Measurement> &dataH,
                                  const bool automated )
{
  const double lowres_min_nsigma_peak = automated ? 2.75 : 1.5;
  const double lowres_min_core_multipeak_chi2dof_peak_improvment = automated ? 0.75 : 0.5;
  const double max_ratio_make_chi2_worse = 2.0;
  const double min_line_chi2 = 1.2;
  const double twopeak_chi2dof_improvement_over_onepeak = 0.5;
  
  if( originalpeaks.empty() )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_multi_peak_fit: should have _some_ originalpeaks!"
         << "  Ignoring this error and continuing" << "\n";
#endif
    
    return RejectedCanTryAgain;
  }
  
  vector<std::shared_ptr<const PeakDef> > toadd = fitpeaks;
  vector<std::shared_ptr<const PeakDef> > toremove = originalpeaks;
    
  while( !toremove.empty() )
  {
    std::shared_ptr<const PeakDef> oldpeak = toremove[0];
    const double mean = oldpeak->mean();
    size_t maxind = 0;
    for( size_t i = 1; i < toadd.size(); ++i )
      if( fabs(toadd[i]->mean()-mean) < fabs(toadd[maxind]->mean()-mean) )
        maxind = i;
    toremove.erase( toremove.begin() );
    toadd.erase( toadd.begin() + maxind );
  }//while( !toremove.empty() )
  
  if( toadd.size() != 1 )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_multi_peak_fit: found toadd.size()=" << toadd.size()
    << " when it should have been 1.  This is a serious logical error,"
    << " but I'm letting it slip and and instead failing these peaks" << "\n";
#endif
    
    return RejectedCanTryAgain;
  }//if( toadd.size() != 1 )
  
  std::shared_ptr<const PeakDef> newpeak = toadd[0];
  const double mean = newpeak->mean();
  const double sigma = newpeak->sigma();
  const double fwhm = newpeak->fwhm();
  const double gausarea = newpeak->gauss_integral( mean-sigma, mean+sigma );
  const double dataarea = dataH->gamma_integral( mean-sigma, mean+sigma );
    
//  cout << "check_lowres_multi_peak_fit: New Peak at mean=" << mean << " and amplitude " << gausarea << endl;
  
  const double newchi2dof = newpeak->chi2dof();
  double oldchi2dof = 0.0;
  for( size_t i = 0; i < originalpeaks.size(); ++i )
    oldchi2dof += originalpeaks[i]->chi2dof();
  oldchi2dof /= originalpeaks.size();  //should actually do based on unique continuums
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  DebugLog(cout) << "New chi2=" << newchi2dof << ", oldchi2=" << oldchi2dof << "\n";
#endif
  
  if( newchi2dof/oldchi2dof > max_ratio_make_chi2_worse )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_multi_peak_fit: Failed add a peak to a ROI for"
    << " for making the chi2DOF much worse "
    << "(newchi2dof:" << newchi2dof << ", oldchi2dof:" << oldchi2dof
    << ")" << "\n";
#endif
    
    return RejectedCanTryAgain;
  }
  
  const double nsigma = gausarea / sqrt( std::max(dataarea,1.0) );
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  DebugLog(cout) << "mean=" << mean << ", nsigma=" << nsigma << "\n";
#endif
  
  if( nsigma < lowres_min_nsigma_peak )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_multi_peak_fit: Failed add a peak to a ROI for"
         << " for not being statistically significant enough "
         << "(nsigma:" << nsigma << ", needed " << lowres_min_nsigma_peak
         << ")" << "\n";
#endif
    
    return RejectedCanTryAgain;
  }
      
  const double core_start = std::max( mean - fwhm, newpeak->lowerX() );
  const double core_end = std::min( mean + fwhm, newpeak->upperX() );
          
  vector<PeakDef> oldpeaks, allpeaks;
  for( size_t i = 0; i < fitpeaks.size(); ++i )
  {
    allpeaks.push_back( *fitpeaks[i] );
    if( fitpeaks[i] != newpeak )
        oldpeaks.push_back( *fitpeaks[i] );
  }//for( size_t i = 0; i < answer.first.size(); ++i )
          
  const double core_chi2dof = evaluate_chi2dof_for_range( allpeaks,
                                                  dataH, core_start, core_end );
  const double no_peak_chi2dof = evaluate_chi2dof_for_range( oldpeaks,
                                                  dataH, core_start, core_end );
          
          
  if( (no_peak_chi2dof - core_chi2dof) < lowres_min_core_multipeak_chi2dof_peak_improvment )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_lowres_multi_peak_fit: Failed add a peak to a ROI for"
         << " for not imroving chi2 enough in center part of new peak "
         << "(noNewPeakChi2Dof:" << no_peak_chi2dof
         << ", withNewPeakChi2Dof:" << core_chi2dof << ")" << "\n";
#endif
    
    return RejectedCanTryAgain;
  }
  
  
  {//begin codeblock to make sure the continuum doesnt dip down all crazy like
    const size_t startchannel = dataH->find_gamma_channel( newpeak->lowerX() );
    size_t endchannel = dataH->find_gamma_channel( newpeak->upperX() );
    const vector<float> &energies = *dataH->gamma_channel_energies();
  
    if( endchannel >= (energies.size()-1) )
      --endchannel;

    size_t minchanel = 0;
    double minval = DBL_MAX;
    std::shared_ptr<const PeakContinuum> continuum = newpeak->continuum();
  
    for( size_t i = startchannel; i <= endchannel; ++i )
    {
      const double val = continuum->offset_integral( energies[i], energies[i+1], dataH );
      if( val < minval )
      {
        minval = val;
        minchanel = i;
      }
    }//for( size_t i = startchannel; i <= endchannel; ++i )
  
    const double lowedgeval = continuum->offset_integral( energies[startchannel],
                                                    energies[startchannel+1], dataH );
    const double highedgeval = continuum->offset_integral( energies[endchannel],
                                                    energies[endchannel+1], dataH );
  
    //THe below 0.5 and 10.0 are based off nearly nothing
    if( minval < 0.5*lowedgeval && minval < 0.5*highedgeval
        && (lowedgeval > 10.0 && highedgeval>10.0) )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_multi_peak_fit: Failed add a peak to a ROI because"
      << " doing so made continuum act badly "
      << "(minval:" << minval << " at channel " << minchanel
      << ", leftval: " << lowedgeval << " at channel " << startchannel
      << " energy " << energies[startchannel]
      << ", rightval: " << highedgeval << " at channel " << endchannel
      << " energy " << energies[endchannel]
      << ")" << "\n";
#endif
      
      return RejectedCanTryAgain;
    }
  }//end codeblock to make sure the continuum doesnt dip down all crazy like
  
  
  //check that we couldnt just draw a line and have the same chi2 as with the
  //  peak
  try
  {
    const vector<float> &energies = *dataH->gamma_channel_energies();
    const vector<float> &contents = *dataH->gamma_channel_contents();
    
    
    const vector<std::shared_ptr<const PeakDef> >::const_iterator pos
                           = find( fitpeaks.begin(), fitpeaks.end(), newpeak );
    if( pos == fitpeaks.end() )
      throw logic_error( "Couldnt find newpeak in fitpeaks" );
    
    const size_t index = pos - fitpeaks.begin();
    double new_peak_lowerx, new_peak_upperx;
    if( index == 0 )
      new_peak_lowerx = newpeak->lowerX();
    else
      new_peak_lowerx = 0.5*(newpeak->mean() + fitpeaks[index-1]->mean());
    
    if( index == (fitpeaks.size()-1) )
      new_peak_upperx = newpeak->upperX();
    else
      new_peak_upperx = 0.5*(newpeak->mean() + fitpeaks[index+1]->mean());
      
    const size_t lowerchannel = dataH->find_gamma_channel( new_peak_lowerx );
    const size_t upperchannel = dataH->find_gamma_channel( new_peak_upperx );
    
    std::vector<double> poly_coeffs, coeff_uncerts;
    const size_t nchannelsrange = (1 + upperchannel - lowerchannel);
    fit_to_polynomial( &(energies[0]) + lowerchannel,
                      &(contents[0]) + lowerchannel,
                      nchannelsrange, 1, poly_coeffs, coeff_uncerts );
    
    const size_t lowercore = dataH->find_gamma_channel( core_start );
    const size_t uppercore = dataH->find_gamma_channel( core_end );
    
    double linechi2 = 0;
    for( size_t channel = lowercore; channel <= uppercore; ++channel )
    {
      const double x = dataH->gamma_channel_lower( channel );
      const double y = dataH->gamma_channel_content( channel );
      const double y_pred = poly_coeffs[0] + poly_coeffs[1]*x;
      const double uncert = (y > MIN_CHANNEL_UNCERT ? sqrt(y) : 1.0);
      linechi2 += std::pow( (y_pred - y) / uncert, 2.0 );
    }//for( int bin = 0; bin < nbin; ++bin )
    
    const double linechi2dof = linechi2 / nchannelsrange;
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "multipeak For " << mean << " core_chi2dof=" << core_chi2dof
    << ", linechi2dof=" << linechi2dof
    << "\n";
#endif
    
    if( linechi2dof < min_line_chi2 )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_multi_peak_fit: Failed add a peak to a ROI for"
      << " for region of peak being compatible with a line "
      << "(linechi2dof:" << linechi2dof
      << ", core_chi2dof:" << core_chi2dof
      << ", new_peak_lowerx: " << new_peak_lowerx
      << ", new_peak_upperx: " << new_peak_upperx
      << ")" << "\n";
#endif
      
      return RejectedCanTryAgain;
    }
    
  }catch( std::exception &e )
  {
    cerr << "check_lowres_multi_peak_fit caught: " << e.what() << endl;
  }
  
  
  //Try to refit ROI using all peaks besides the one fit, to see if the chi2dof
  //  improves with having the additional peak, verses going form a linear to
  //  quadratic continuum,

//  if( fitpeaks.size() == 3 )
  {
    vector< std::shared_ptr<const PeakDef> > otherpeak;
    for( size_t i = 0; i < fitpeaks.size(); ++i )
      if( fitpeaks[i] != newpeak )
        otherpeak.push_back( fitpeaks[i] );
  
    const double lx = newpeak->lowerX();
    const double ux = newpeak->upperX();
    double withoutNewChi2Dof;
    vector<PeakDef> withoutResultPeaks;
    const double withNewChi2Dof = newpeak->chi2dof();
    
    refit_for_new_roi( otherpeak, dataH, lx, ux, withoutNewChi2Dof, withoutResultPeaks );
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "second peak at " << newpeak->mean() << " original peak at "
         << otherpeak[0]->mean() << ", withNewChi2Dof=" << withNewChi2Dof
         << ", withoutNewChi2Dof=" << withoutNewChi2Dof
         << ", diff=" << (withoutNewChi2Dof - withNewChi2Dof) << "\n";
#endif
    
    if( withoutResultPeaks.size()
        && (withoutNewChi2Dof - withNewChi2Dof) < twopeak_chi2dof_improvement_over_onepeak )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_lowres_multi_peak_fit: Failed add a peak at "
           << newpeak->mean() << " to a ROI (for peak at "
           << otherpeak[0]->mean() << ") for not giving a big enough"
           << " improvement in chi2, vs just incrementing the order of"
           << " polynomial continuum." << "\n";
#endif
      
      return RejectedDontTryAgain;
    }
  }//if( fitpeaks.size() == 2 )
  
  
  return AcceptPeak;
}//PeakRejectionStatus check_lowres_multi_peak_fit(...)



PeakRejectionStatus check_highres_multi_peak_fit( const vector<std::shared_ptr<const PeakDef> > fitpeaks,
                                                const vector<std::shared_ptr<const PeakDef> > originalpeaks,
                                                const std::shared_ptr<const Measurement> &dataH,
                                                const bool automated )
{
  typedef std::shared_ptr<const PeakDef> PeakDefShrdPtr;
  
  const double min_sigma_nearest_existing_peak = 1.0;
  const double max_chi2_dof = 25.0;
  const double max_nsigma_newpeak_to_nearest_other = 10.0; //guessed
  const double min_nsigma_peak = automated ? 3.0 : 1.5;
  
  vector<std::shared_ptr<const PeakDef> > toadd = fitpeaks;
  vector<std::shared_ptr<const PeakDef> > toremove = originalpeaks;

  
  while( !toremove.empty() )
  {
    std::shared_ptr<const PeakDef> oldpeak = toremove[0];
    const double mean = oldpeak->mean();
    size_t maxind = 0;
    for( size_t i = 1; i < toadd.size(); ++i )
      if( fabs(toadd[i]->mean()-mean) < fabs(toadd[maxind]->mean()-mean) )
        maxind = i;
    toremove.erase( toremove.begin() );
    toadd.erase( toadd.begin() + maxind );
  }//while( !toremove.empty() )
  
  if( toadd.size() != 1 )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_highres_multi_peak_fit: found toadd.size()=" << toadd.size()
    << " when it should have been 1.  This is a serious logical error,"
    << " but I'm letting it slip and and instead failing these peaks" << "\n";
#endif
    
    return RejectedCanTryAgain;
  }//if( toadd.size() != 1 )
  
  std::shared_ptr<const PeakDef> newpeak = toadd[0];
  const double mean = newpeak->mean();
  const double sigma = newpeak->sigma();
//  const double fwhm = newpeak->fwhm();
  const double gausarea = newpeak->gauss_integral( mean-sigma, mean+sigma );
  const double dataarea = dataH->gamma_integral( mean-sigma, mean+sigma );
  
  const bool debug_this_peak = false; //fabs(mean - 729.293) < 2.0;
  
  if( debug_this_peak )
    cout << "debug_this_peak" << endl;
  
  
  const double newchi2dof = newpeak->chi2dof();
  double oldchi2dof = 0.0;
  for( size_t i = 0; i < originalpeaks.size(); ++i )
    oldchi2dof += originalpeaks[i]->chi2dof();
  oldchi2dof /= originalpeaks.size();  //should actually do based on unique continuums
  
  //Check that the ROI Chi2 isnt horrible
  if( newchi2dof > max_chi2_dof && newchi2dof > oldchi2dof )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_highres_multi_peak_fit: Failed to fit a peak for it having "
         << " too high a chi2dof (" << newchi2dof << ", was " << oldchi2dof
         << ")" << "\n";
#endif
    return RejectedCanTryAgain;
  }//if( newchi2dof > max_chi2_dof && newchi2dof > oldchi2dof )
  
  //Check that the new peak isnt super close to another peak
  for( const PeakDefShrdPtr &p : fitpeaks )
  {
    if( p == newpeak )
      continue;
    
    const double dx = fabs( p->mean() - mean );
    
    if( (dx/sigma) < min_sigma_nearest_existing_peak )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_highres_multi_peak_fit: Failed to fit a peak for it being "
           << " to close to another peak (" << mean << " vs " << p->mean()
           << " with sigma=" << sigma << "})" << "\n";
#endif
      return RejectedDontTryAgain;
    }
  }//for( const PeakDefShrdPtr &p : fitpeaks )
  
  
  //Check that none of the peaks have gone outside of the sigma range they
  //  should stay in.
  for( const PeakDefShrdPtr &p : fitpeaks )
  {
    const float mean = static_cast<float>( p->mean() );
    const float sigma = static_cast<float>( p->sigma() );
    
    float min_sigma, max_sigma;
    expected_peak_width_limits( mean, true, dataH, min_sigma, max_sigma );
    
    bool outsideExpectedFwhm = (sigma < min_sigma || sigma > max_sigma);
    
    // We checked against reasonable expected FWHM, but incase this failed for some reason,
    //  we'll give it another opportunity by seeing if channel counts are reasonable.
    if( outsideExpectedFwhm )
    {
      auto cal = dataH->energy_calibration();
      if( cal && cal->valid() )
      {
        // TODO: The valid number of channels in a peak has not been looked at very closely at all
        const double min_num_channel = 2.5;  // 2.5 seen on HPGe with 4096 channels and 3 MeV scale
        const double max_num_channel = 15.0; // 12 seen for the 2614 of a 16k channel, 3 MeV spec
        
        const double lowerSigmaChannel = cal->channel_for_energy( mean - sigma );
        const double upperSigmaChannel = cal->channel_for_energy( mean + sigma );
        const double nchandiff = upperSigmaChannel - lowerSigmaChannel;
        outsideExpectedFwhm = (nchandiff > min_num_channel && nchandiff < max_num_channel);
      }
    
    }//if( outsideExpectedFwhm )
    
    if( outsideExpectedFwhm )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_highres_multi_peak_fit: Failed to fit a peak for the width"
      << " being out of range (" << sigma << ", expected {"
      << min_sigma << ", " << max_sigma << "})" << "\n";
#endif
      
      return RejectedDontTryAgain;
    }
  }//for( const PeakDefShrdPtr &p : fitpeaks )
  
  
  
  //reject if a peak is more than X sigma from nearest peak.
  double smallestdx = DBL_MAX;
  for( const PeakDefShrdPtr &p : fitpeaks )
  {
    if( p != newpeak )
      smallestdx = std::min( fabs(p->mean() - newpeak->mean()), smallestdx );
  }//for( const PeakDefShrdPtr &p : fitpeaks )
  
//  if( mean < 100.0 )
//    cerr << "For mean " << mean << " (smallestdx/sigma)=" << (smallestdx/sigma) << endl;
  
  if( (smallestdx/sigma) > max_nsigma_newpeak_to_nearest_other )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_highres_multi_peak_fit: Failed to fit a peak for it being "
         << " to far away from nearest neighbor (dist=" << smallestdx
         << ", sigma=" << sigma << ")" << "\n";
#endif
    return RejectedCanTryAgain;
  }//
  
  
  //reject if any mean is outside ROI  (RejectedCanTryAgain)
  for( const PeakDefShrdPtr &p : fitpeaks )
  {
    if( p->mean() < p->lowerX() || p->mean() > p->upperX() )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cerr) << "check_highres_multi_peak_fit: Failed to fit a peak for causing "
           << " a peak sharing the continuum to have a mean outside the ROI "
           << "(mean=" << p->mean() << ", roilow=" << p->lowerX()
           << ", roiupper=" << p->upperX() << ")" << "\n";
#endif
      return RejectedCanTryAgain;
    }
  }//for( const PeakDefShrdPtr &p : fitpeaks )
  
  
  //reject if a peak is not statistically significant (RejectedCanTryAgain)
  const double nsigma = gausarea / sqrt( std::max(dataarea,1.0) );
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  DebugLog(cout) << "mean=" << mean << ", nsigma=" << nsigma << "\n";
#endif
  
  if( nsigma < min_nsigma_peak )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cerr) << "check_highres_multi_peak_fit: Failed add a peak to a ROI for"
         << " for not being statistically significant enough "
         << "(nsigma:" << nsigma << ", needed " << min_nsigma_peak
         << ")" << "\n";
#endif
    
    return RejectedCanTryAgain;
  }//if( nsigma < min_nsigma_peak )
  
  return AcceptPeak;
}//PeakRejectionStatus check_highres_multi_peak_fit(...)



bool check_highres_single_peak_fit( const std::shared_ptr<const PeakDef> peak,
                                   const std::shared_ptr<const Measurement> &dataH,
                                   const bool automated )
{
  const double max_chi2dof_persignificance = 2.0;
  const double min_nsigma_peak = automated ? 5.0 : 1.5;
  const double low_stat_min_nsigma_peak = automated ? 3.5 : 0.5;
  const double med_stat_min_nsigma_peak = automated ? 4.25 : 0.75;
  const double min_core_chi2dof_peak_improvment = automated ? 0.5 : 0.25;
  const double min_significance_test_width = 7.5;
  const double min_chi2dof_test_width = 3.0;
  const double max_chi2dof_roi = automated ? 50.0 : 250.0;
  const double max_nsignif_to_apply_chidof_check = 100.0; //based on a single example (80.0 keV peak of example Ba133 spectrum)
  
  const double fwhm = peak->fwhm();
  const double mean = peak->mean();
  const double sigma = peak->sigma();
  const double chi2dof = peak->chi2dof();
  const double core_start = std::max( mean - fwhm, peak->lowerX() );
  const double core_end = std::min( mean + fwhm, peak->upperX() );
  

  const bool debug_this_peak = false; //fabs(mean - 32.9752) < 2.5;
  
  if( debug_this_peak )
    cout << "debug_this_peak" << endl;
  
  vector<PeakDef> fitpeaks( 1, *peak ), fitpeaksnoamp( 1, *peak );
  fitpeaksnoamp[0].setAmplitude( 0.0 );
  const double core_chi2dof = evaluate_chi2dof_for_range( fitpeaks,
                                                         dataH, core_start, core_end );
  const double no_peak_chi2dof = evaluate_chi2dof_for_range( fitpeaksnoamp,
                                                            dataH, core_start, core_end );
  
  // Check to see if having the peak be there is actually any improvement to the fit.
  //  However, if the Chi/dof is already extremely good, just require there is some improvement
  //  (this happens on artificially good statistics spectra)
  // TODO: The "4" below is totally arbitrary - should have some sliding scale or something
  if( (core_chi2dof > (4*min_core_chi2dof_peak_improvment)) )
  {
    if( (no_peak_chi2dof - core_chi2dof) < min_core_chi2dof_peak_improvment )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      if( debug_this_peak )
        DebugLog(cerr) << "check_highres_single_peak_fit: Failed to fit a peak for the "
        << "change in chi2 without the peak being to small (chi2dof="
        << no_peak_chi2dof << " w/ continuum only, " << core_chi2dof
        << " w/ continuum and peak)" << "\n";
#endif
      
      return false;
    }//if( (no_peak_chi2dof - core_chi2dof) < 0.5 )
  }else
  {
    if( (no_peak_chi2dof - core_chi2dof) < 0.0 )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      if( debug_this_peak )
        DebugLog(cerr) << "check_highres_single_peak_fit: Failed to fit a peak because "
        << "the continuum only chi2 (chi2dof="
        << no_peak_chi2dof << ") was better than with the peak (chi2dof=" << core_chi2dof
        << " w/ continuum and peak)\n";
#endif
      
      return false;
    }//if( (no_peak_chi2dof - core_chi2dof) < 0.5 )

  }//if( (core_chi2dof > (4*min_core_chi2dof_peak_improvment)) )
  
  
  const double gausarea = peak->gauss_integral( mean-2.0*sigma, mean+2.0*sigma );
  const double dataarea = dataH->gamma_integral( mean-2.0*sigma, mean+2.0*sigma );
  const double nsigma = gausarea / sqrt( std::max(dataarea,1.0) );
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  if( debug_this_peak )
  {
    const double chi2Dof = peak->chi2dof();
    DebugLog(cerr) << "mean: " << mean << ", nsigma=" << nsigma  << ", chi2Dof/nsigma="
         << chi2Dof/nsigma << "\n";
  }//
#endif
  
  if( (peak->chi2dof()/nsigma) > max_chi2dof_persignificance )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    if( debug_this_peak )
      DebugLog(cerr) << "check_highres_single_peak_fit: Failed to fit a peak for the chi2DOF"
           << " being to bad (chi2dof=" << peak->chi2dof()
           << ", nsigma_significance=" << nsigma << ")" << "\n";
#endif
    return false;
  }//if( chi2Dof > 25.0 )
  
  const double nsignif = peak->amplitude() / peak->amplitudeUncert();
  const bool lowstatregion = ((dataarea - gausarea) <= 3.0*sqrt(dataarea));
  const bool medstatregion = ((dataarea - gausarea) <= 9.0*sqrt(dataarea));

  if( nsignif < min_nsigma_peak
     && (!lowstatregion || nsignif < low_stat_min_nsigma_peak || peak->amplitude()<15.0)
     && (!medstatregion || nsignif < med_stat_min_nsigma_peak || peak->amplitude()<25.0) )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    if( debug_this_peak )
      DebugLog(cerr) << "check_highres_single_peak_fit: Failed to fit a peak for the gross "
           << "stat significance to low (" << nsigma << ")" << "\n";
#endif
    
    return false;
  }
  
  
  
  if( chi2dof > max_chi2dof_roi && nsignif < max_nsignif_to_apply_chidof_check )
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    if( debug_this_peak )
      DebugLog(cerr) << "check_highres_single_peak_fit: Failed to fit a peak chi2dof "
      << "too high (mean=" << mean << ", chi2dof=" << chi2dof
      << ", nsignif=" << nsignif << ")\n";
#endif
    
    return false;
  }
  
  
  float min_sigma, max_sigma;
  expected_peak_width_limits( (float)mean, true, dataH, min_sigma, max_sigma );
  
  //An issue is that doppler broadened peaks (like 511 keV) will have a width
  //  outside of limits - so if the chi2dof is good enough, or the peak is
  //  really significant, we'll let this test slide
  if( (nsignif < min_significance_test_width || core_chi2dof > min_chi2dof_test_width)
      && (sigma < min_sigma || sigma > max_sigma) )
  {
    if( debug_this_peak || PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      cerr << "check_highres_single_peak_fit: Failed to fit a peak for the width"
           << " being out of range (mean=" << mean << ", sigma=" << sigma
           << ", expected {" << min_sigma << ", " << max_sigma << "})" << "\n";
    
    return false;
  }
  
  
  //if we're near the turn on of the detector, we'll try to
  const size_t meanchannel = dataH->find_gamma_channel( mean );
  const size_t nchannel = dataH->num_gamma_channels();
  const double fracchannel = static_cast<double>(meanchannel) / nchannel;
  const double turnonfraction = 0.03;  //90 keV on 3 MeV scale, 240 keV on 8 MeV scale
  if( automated && fracchannel < turnonfraction )
  {
    size_t lower_channel, upper_channel;
    const bool hasextent = ExperimentalPeakSearch::find_spectroscopic_extent(
                                          dataH, lower_channel, upper_channel );

    const float lowextent = dataH->gamma_channel_lower( lower_channel );
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    cout << "Lower extent=" << lowextent << endl;
#endif
    
    const size_t minus1Sigmachannel = dataH->find_gamma_channel( mean - sigma );
    
    
    if( hasextent && (minus1Sigmachannel < lower_channel) )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      if( debug_this_peak )
        DebugLog(cerr) << "check_highres_single_peak_fit: Failed to fit a because the"
        << " mean is below spectroscopic extent for this spectrum"
        << " (mean=" << mean << ", channel(mean-sigma)=" << minus1Sigmachannel
        << ", lowerspecchannel=" << lower_channel <<  ")" << "\n";
#endif
      
      return false;
    }
    
/*
    const double clx = peak->lowerX();
    const double cux = peak->upperX();
    const size_t lxchannel = dataH->find_gamma_channel( clx );
    const size_t uxchannel = dataH->find_gamma_channel( cux );
    const size_t m2sigmachannel = dataH->find_gamma_channel( mean-2.0*sigma );
  
    const float lxl = dataH->gamma_channel_lower(lxchannel);
    const float lxu = dataH->gamma_channel_upper(lxchannel);
    const float uxl = dataH->gamma_channel_lower(uxchannel);
    const float uxu = dataH->gamma_channel_upper(uxchannel);
    const float mxl = dataH->gamma_channel_lower(meanchannel);
    const float mxu = dataH->gamma_channel_upper(meanchannel);
    const float n2xl = dataH->gamma_channel_lower(m2sigmachannel);
    const float n2xu = dataH->gamma_channel_upper(m2sigmachannel);
    
    
    
    const double lcont = peak->continuum()->offset_integral( lxl, lxu );
    const double ucont = peak->continuum()->offset_integral( uxl, uxu );
    const double cslope = (ucont - lcont) / (uxu - lxl);
    
    const double mp = peak->gauss_integral( mxl, mxu );
    const double n2p = peak->gauss_integral( n2xl, n2xu );
    const double pslope = (mp - n2p) / (n2xu - n2xl);
    
    cout << mean << ", cslope=" << cslope << ", pslope=" << pslope << "-->" << cslope/pslope << endl;

    
    cout << "lower_channel=" << lower_channel << ", upper_channel=" << upper_channel << ", meanchannel=" << meanchannel << endl;
    
    
    const double max_ratio_continuum_to_peak_slope = 1.0;
    if( false && fabs(cslope/pslope) > max_ratio_continuum_to_peak_slope )
    {
      if( debug_this_peak || PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cerr) << "check_highres_single_peak_fit: Failed to fit a because the"
             << " continuum slope is too steep relative to the peak slope"
             << " (mean=" << mean << ", continuumslope=" << cslope
             << ", peakslope(2sig)=" << pslope <<  ")" << "\n";

      return false;
    }//if( fabs(cslope/pslope) > max_ratio_continuum_to_peak_slope )
*/
  }//if( invfracchannel < invturnonfraction )
  
  
  
  //check that we couldnt just draw a line and have the same chi2 as with the
  //  peak
  try
  {
    const vector<float> &energies = *dataH->gamma_channel_energies();
    const vector<float> &contents = *dataH->gamma_channel_contents();
    const size_t lowerchannel = dataH->find_gamma_channel( peak->lowerX() );
    const size_t upperchannel = dataH->find_gamma_channel( peak->upperX() );
    
    std::vector<double> poly_coeffs, coeff_uncerts;
    const size_t nchannelsrange = (1 + upperchannel - lowerchannel);
    fit_to_polynomial( &(energies[0]) + lowerchannel,
                      &(contents[0]) + lowerchannel,
                      nchannelsrange, 1, poly_coeffs, coeff_uncerts );
    
    const size_t lowercore = dataH->find_gamma_channel( core_start );
    const size_t uppercore = dataH->find_gamma_channel( core_end );
    
    double linechi2 = 0.0, avrguncert = 0.0;
    for( size_t channel = lowercore; channel <= uppercore; ++channel )
    {
      const double x = dataH->gamma_channel_lower( channel );
      const double y = dataH->gamma_channel_content( channel );
      const double y_pred = poly_coeffs[0] + poly_coeffs[1]*x;
      const double uncert = (y > MIN_CHANNEL_UNCERT ? std::sqrt(y) : 1.0);
      const double thichi2 = std::pow( (y_pred - y) / uncert, 2.0 );
      
      linechi2 += thichi2;
      avrguncert += 1.0 / uncert;
    }//for( int bin = 0; bin < nbin; ++bin )
    
    avrguncert /= (uppercore - lowercore + 1);
    const double linechi2dof = linechi2 / (uppercore - lowercore + 1);
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    if( debug_this_peak )
      DebugLog(cout) << "For " << mean << " core_chi2dof=" << core_chi2dof
           << ", linechi2dof=" << linechi2dof << ", nchannelsrange="
           << nchannelsrange << "\n"
           << "avrguncert=" << avrguncert << "\n";
#endif
    
    const double n1sigma = mean - sigma;
    const double p1sigma = mean + sigma;
    const double peakarea = peak->gauss_integral( n1sigma, p1sigma );
    const double dataarea = dataH->gamma_integral( n1sigma, p1sigma );
    
    size_t lbin = dataH->find_gamma_channel( n1sigma );
    size_t rbin = dataH->find_gamma_channel( p1sigma );
    const size_t nbin_pm1sigma = (rbin - lbin + 1);
    
    const double peak_above_uncert = peakarea / sqrt( dataarea );
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    if( debug_this_peak )
      DebugLog(cout) << "peak_above_uncert=" << peak_above_uncert
           << " over nbins=" << nbin_pm1sigma
           << "--->" << (peak_above_uncert/nbin_pm1sigma) << "\n"
           << "For " << mean
           << "\tcore_chi2dof=" << core_chi2dof << "\n"
           << "\tlinechi2dof=" << linechi2dof << "\n"
           << "\tno_peak_chi2dof=" << no_peak_chi2dof << "\n"
           << "\t(linechi2dof-core_chi2dof)=" << (linechi2dof - core_chi2dof) << "\n"
           << "\t(linechi2dof-core_chi2dof)/linechi2dof=" << (linechi2dof - core_chi2dof)/linechi2dof << "\n"
           << "\t(core_chi2dof/linechi2dof)=" << (core_chi2dof/linechi2dof) << "\n"
           << "\t(linechi2dof/core_chi2dof)=" << (linechi2dof/core_chi2dof)
           << ", (no_peak_chi2dof/linechi2dof)=" << (no_peak_chi2dof/linechi2dof)
           << "\n";
#endif
    
    //None of the below cuts have been validated for highres spectra
    const double min_core_chi2dof_over_line_improvment  = automated ? 0.8 : 0.5;
    const double min_core_line_chi2dof_lower_ratio = automated ? 0.5 : 0.25;
    const double bad_continuum_fit_multiple = automated ? 0.333 : 0.0;
    const double min_necassary_chi2_improv_over_line = automated ? 0.35 : 0.1;
    const double max_avrg_uncert_to_require_chi2_improv_over_line = 0.45;
    const double badcont_line_chi2dof = 4.0;
    const double badcont_corecontinuum_chi2dof = 10.0;
    const double badcont_cont_to_line_ratio = 3.5;
    const double line_is_close = 1.1;
    const double line_is_close_necassary_improvment = 0.1;
    const double badcont_upper_chi2diff = 0.75;
    const double max_line_chi2dof = automated ? 1.5 : 0.25;
    const double medstat_max_line_chi2dof = automated ? 1.25 : 0.25;
    
    if( linechi2dof < max_line_chi2dof
        && (!medstatregion || linechi2dof < medstat_max_line_chi2dof)
        && !lowstatregion )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      if( debug_this_peak )
        DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
             << "data was consistent with a straight line"
             << " (linechi2dof=" << linechi2dof << ")" << "\n";
#endif
      return false;
    }
    
    //if the continuum only chi2 is much different then the line-only chi2,
    if( (linechi2dof - core_chi2dof) < min_core_chi2dof_over_line_improvment
       && ((core_chi2dof/linechi2dof) > min_core_line_chi2dof_lower_ratio)
       && ((peak_above_uncert/nbin_pm1sigma) > 0.85
           || (linechi2dof/core_chi2dof) < bad_continuum_fit_multiple*(no_peak_chi2dof/linechi2dof)) )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      if( debug_this_peak )
        DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
             << "change in chi2 over a straight line being to small (core_chi2dof="
             << core_chi2dof << ", linechi2dof=" << linechi2dof << ")"
             << "\n";
#endif
      
      return false;
    }//if( (no_peak_chi2dof - core_chi2dof) < 0.5 )
    
    if( (linechi2dof-core_chi2dof) < min_necassary_chi2_improv_over_line
       && avrguncert < max_avrg_uncert_to_require_chi2_improv_over_line )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      if( debug_this_peak )
        DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
             << "change in chi2 over a straight line being to small no matter what"
             << " (core_chi2dof=" << core_chi2dof << ", linechi2dof="
             << linechi2dof << ")" << "\n";
#endif
      return false;
    }
    
    //This next test is intended to elimnate peaks where the continuum is
    //  much below where it visually should be, but because of the gausians
    //  amplitude the overall chi2 is decent in the fit with the peak, but
    //  really this is due to having the incorrect continuum.
    if( linechi2dof < badcont_line_chi2dof
       && no_peak_chi2dof > badcont_corecontinuum_chi2dof
       && (no_peak_chi2dof/linechi2dof) > badcont_cont_to_line_ratio
       && (linechi2dof-core_chi2dof) < badcont_upper_chi2diff )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      if( debug_this_peak )
        DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
             << "a line fit the region well, but the peaks continuum didnt"
             << " (core_chi2dof=" << core_chi2dof  << ", no_peak_chi2dof="
             << no_peak_chi2dof << ", linechi2dof=" << linechi2dof
             << ")" << "\n";
#endif
      return false;
    }
    
    if( linechi2dof < line_is_close
       && (linechi2dof-core_chi2dof) < line_is_close_necassary_improvment )
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      if( debug_this_peak )
        DebugLog(cerr) << "check_lowres_single_peak_fit: Failed to fit a peak because the "
             << "a line fit the region well, and adding a peak didnt improve fit enough"
             << " (linechi2dof=" << linechi2dof  << ", core_chi2dof="
             << core_chi2dof << ")" << "\n";
#endif
      return false;
    }
  }catch( std::exception &e )
  {
    cerr << "check_highres_single_peak_fit(): caught " << e.what() << endl;
  }//try / catch
  
  
  return true;
}//bool check_highres_single_peak_fit(...)





pair< PeakShrdVec, PeakShrdVec > searchForPeakFromUser( const double x,
                                                        double pixelPerKev,
                                                        const std::shared_ptr<const Measurement> &dataH,
                                                        const PeakShrdVec &inpeaks,
                                                        std::shared_ptr<const DetectorPeakResponse> drf,
                                                       const bool isHPGe )
{
  typedef std::shared_ptr<const PeakDef> PeakDefShrdPtr;
  
  if( !dataH || !dataH->num_gamma_channels() )
    return pair<PeakShrdVec,PeakShrdVec>();
  
  const bool automated = (pixelPerKev <= 0.0);
  if( automated )
    pixelPerKev = 5.0;
  
  const bool lowres_shrink_roi_for_multipeaks = true;
  
  bool lowstatregion = false;
  const size_t nchannels = dataH->num_gamma_channels();
  
  double sigma0, mean0, area0;
  get_candidate_peak_estimates_for_user_click( sigma0, mean0, area0, x,
                                              pixelPerKev, dataH, isHPGe, inpeaks );
  
  if( drf && drf->isValid() && drf->hasResolutionInfo() )
    sigma0 = drf->peakResolutionSigma( mean0 );
  
  PeakShrdVec coFitPeaks;
  double roiLower, roiUpper;
  
  try
  {
    combine_peaks_to_roi( coFitPeaks, roiLower, roiUpper, lowstatregion,
                          dataH, inpeaks, mean0, sigma0, area0, pixelPerKev, isHPGe );
  }catch( std::exception &e )
  {
    if( !inpeaks.empty() )
    {
      return pair<PeakShrdVec,PeakShrdVec>();
    }else
    {
      const double rough_fwhm = isHPGe ? PeakFitUtils::hpge_fwhm_fcn( static_cast<float>(mean0) )
                                        : PeakFitUtils::nai_fwhm_fcn( static_cast<float>(mean0) );
      cerr << "Failed to calc ROI range, will WAG" << endl;
      roiLower = mean0 - std::max( 2*std::max(rough_fwhm,sigma0), 3.0 );
      roiUpper = mean0 + std::max( 2*std::max(rough_fwhm,sigma0), 3.0 );
    }
  }// try / catch
  
  const size_t nFitPeaks = coFitPeaks.size() + 1;
  
  if( coFitPeaks.size() )
  {
    //we should update sigma 0 to be for lowest energy peak, and start out a
    //  bit narrower (this is purely a guess at this point, we should really
    //  do something that takes into account sigma relative to the ROI or
    //  something)
    if( coFitPeaks.size()==1 && fabs(coFitPeaks[0]->mean()-x) < 2.0*sigma0 )
      sigma0 = coFitPeaks[0]->sigma() / 2.0;
    else
      sigma0 = coFitPeaks[0]->sigma() * nFitPeaks / (nFitPeaks + 1);
    mean0 = x;
  }//if( coFitPeaks.size() > 1 )
  
  
  double chi2Dof;
  vector<double> lowerEnergies, upperEnergies;
  lowerEnergies.push_back( roiLower );
  upperEnergies.push_back( roiUpper );
  
  /** Using google Ceres to fit peaks is an experiment in using Ceres.
   
   Its actually maybe a slower (although this is probably not the fault of Ceres - but of my
   coding - and actually this is only when using a single thread, and for a single peak - for
   multiple peaks and multiple Ceres threads, it doe beat-out minuit - even with my non-optimal
   usage of Ceres), is totally untested, and hasnt been finished implementing everything.
   
   */
#if( USE_REL_ACT_TOOL )
#ifdef _MSC_VER
#pragma message( "Not using L-M peak fit, even though USE_REL_ACT_TOOL defined." )
#else
#warning "Not using L-M peak fit, even though USE_REL_ACT_TOOL defined."
#endif
#define USE_LM_PEAK_FIT 0
#else
#define USE_LM_PEAK_FIT 0
#endif

  
#if( !USE_LM_PEAK_FIT )
  PeakShrdVec initialfitpeaks;
  fit_peak_for_user_click( initialfitpeaks, chi2Dof, dataH, coFitPeaks,
                          mean0, sigma0, area0, lowerEnergies, upperEnergies, isHPGe );
#else
  PeakShrdVec mnInitialfitpeaks;
  const auto t1 = std::chrono::high_resolution_clock::now();
  fit_peak_for_user_click( mnInitialfitpeaks, chi2Dof, dataH, coFitPeaks,
                           mean0, sigma0, area0, lowerEnergies, upperEnergies, isHPGe );
  const auto t2 = std::chrono::high_resolution_clock::now();
  
  for( size_t i = 0; i < mnInitialfitpeaks.size(); ++i )
  {
    cout << "PRE LM  Peak " << std::setw(2) << i << ": mean=" << std::setw(10) << mnInitialfitpeaks[i]->mean()
    << " keV, FWHM=" << std::setw(10) << mnInitialfitpeaks[i]->fwhm() << ", amp=" << std::setw(10)
    << mnInitialfitpeaks[i]->amplitude() << endl;
  }
  
  PeakShrdVec lmInitialfitpeaks;
  const auto t3 = std::chrono::high_resolution_clock::now();
  PeakFitLM::fit_peak_for_user_click_LM( lmInitialfitpeaks, chi2Dof, dataH, coFitPeaks,
                             mean0, sigma0, area0, lowerEnergies[0], upperEnergies[0], isHPGe );
  const auto t4 = std::chrono::high_resolution_clock::now();
  
  
  cout << "Old way fit " << mnInitialfitpeaks.size() << " peaks - LM way fit " << lmInitialfitpeaks.size() << endl;
  const size_t npeaks = std::min( lmInitialfitpeaks.size(), mnInitialfitpeaks.size() );
  
  auto mnCont = mnInitialfitpeaks.size() ? mnInitialfitpeaks[0]->continuum() : nullptr;
  auto lmCont = lmInitialfitpeaks.size() ? lmInitialfitpeaks[0]->continuum() : nullptr;
  
  cout << "Continuum: minuit={";
  if( mnCont )
  {
    for( size_t i = 0; i < mnCont->parameters().size(); ++i )
      cout << (i ? ", " : "") << mnCont->parameters()[i] << " +- " << (i < mnCont->uncertainties().size() ? mnCont->uncertainties()[i] : 0.0);
  }
  cout << "}, L-M={";
  if( lmCont )
  {
    for( size_t i = 0; i < lmCont->parameters().size(); ++i )
    cout << (i ? ", " : "") << lmCont->parameters()[i] << " +- " << (i < lmCont->uncertainties().size() ? lmCont->uncertainties()[i] : 0.0);
  }
  cout << "}\n";
  
  
  for( size_t i = 0; i < npeaks; ++i )
  {
    const PeakDef &mnPeak = *mnInitialfitpeaks[i];
    const PeakDef &lmPeak = *lmInitialfitpeaks[i];
    cout << "Peak " << i << ":\n";
    cout << "\tMinuit" << std::setw(24) << "" << "\t\tL-M" << endl;
    cout << "\tmean: " << std::setw(10) << mnPeak.mean() << " +- " << std::setw(10) << mnPeak.meanUncert()
    << "\t\t" << std::setw(10) << lmPeak.mean() << " +- " << std::setw(10) << lmPeak.meanUncert() << endl;
    cout << "\tsigma: " << std::setw(10) << mnPeak.sigma() << " +- " << std::setw(10) << mnPeak.sigmaUncert()
    << "\t\t" << std::setw(10) << lmPeak.sigma() << " +- " << std::setw(10) << lmPeak.sigmaUncert() << endl;
    cout << "\tamp:   " << std::setw(10) << mnPeak.amplitude() << " +- " << std::setw(10) << mnPeak.amplitudeUncert()
    << "\t\t" << std::setw(10) << lmPeak.amplitude() << " +- " << std::setw(10) << lmPeak.amplitudeUncert() << endl;
    cout << "\tchi2:  " << std::setw(24) << mnPeak.chi2dof() << "\t\t" << std::setw(24) << lmPeak.chi2dof() << endl;
  }
  
  cout << "Old method: " << std::chrono::duration<double, std::milli>(t2 - t1).count()
  << "ms, vs new method " << std::chrono::duration<double, std::milli>(t4 - t3).count()
  << endl;
  
  
  PeakShrdVec initialfitpeaks = lmInitialfitpeaks;
#endif // !USE_LM_PEAK_FIT / else
  
  if( initialfitpeaks.empty() )
    return pair<PeakShrdVec,PeakShrdVec>();
  
  if( initialfitpeaks.size() > 1 )
  {
    PeakRejectionStatus status;
    
    if( isHPGe )
      status = check_highres_multi_peak_fit( initialfitpeaks, coFitPeaks,
                                                 dataH, automated );
    else
      status = check_lowres_multi_peak_fit( initialfitpeaks, coFitPeaks,
                                           dataH, automated );

    switch( status )
    {
      case AcceptPeak:
        return make_pair( initialfitpeaks, coFitPeaks );
        
      case RejectedCanTryAgain:
        break;
        
      case RejectedDontTryAgain:
        return pair<PeakShrdVec,PeakShrdVec>();
    }//switch( status )
      
    
    const PeakDef &lpeak = *coFitPeaks.front();
    const PeakDef &rpeak = *coFitPeaks.back();
    
    if( mean0 < (lpeak.mean()-3.0*lpeak.sigma())
        || mean0 > (rpeak.mean()+3.0*rpeak.sigma()))
    {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
      DebugLog(cout) << "mean0=" << mean0 << ", (mean0-lpeak.mean())/lpeak.sigma()="
           << (mean0-lpeak.mean())/lpeak.sigma() << "\n"
           << "mean0=" << mean0 << ", (mean0-rpeak.mean())/rpeak.sigma()="
           << (mean0-rpeak.mean())/rpeak.sigma() << "\n";
#endif
      
      try
      {
        combine_peaks_to_roi( coFitPeaks, roiLower, roiUpper, lowstatregion,
                    dataH, PeakShrdVec(), mean0, sigma0, area0, pixelPerKev, isHPGe );
      }catch( std::exception & )
      {
        return pair<PeakShrdVec,PeakShrdVec>();
      }
      
      lowerEnergies = vector<double>( 1, roiLower );
      upperEnergies = vector<double>( 1, roiUpper );
      
      fit_peak_for_user_click( initialfitpeaks, chi2Dof, dataH, coFitPeaks,
                          mean0, sigma0, area0, lowerEnergies, upperEnergies, isHPGe );
      
#if( USE_LM_PEAK_FIT )
      for( size_t i = 0; i < initialfitpeaks.size(); ++i )
      {
        cout << "OLD Peak " << std::setw(2) << i << ": mean=" << std::setw(10) << initialfitpeaks[i]->mean()
        << " keV, FWHM=" << std::setw(10) << initialfitpeaks[i]->fwhm() << ", amp=" << std::setw(10)
        << initialfitpeaks[i]->amplitude() << endl;
      }
      
      
      PeakFitLM::fit_peak_for_user_click_LM( initialfitpeaks, chi2Dof, dataH, coFitPeaks,
                                 mean0, sigma0, area0, lowerEnergies[0], upperEnergies[0], isHPGe );
      
      for( size_t i = 0; i < initialfitpeaks.size(); ++i )
      {
        cout << "LM  Peak " << std::setw(2) << i << ": mean=" << std::setw(10) << initialfitpeaks[i]->mean()
        << " keV, FWHM=" << std::setw(10) << initialfitpeaks[i]->fwhm() << ", amp=" << std::setw(10)
        << initialfitpeaks[i]->amplitude() << endl;
      }
#endif  //!USE_LM_PEAK_FIT / else
    }else
    {
      return pair<PeakShrdVec,PeakShrdVec>();
    }
  }//if( !isHPGe && initialfitpeaks.size() > 1 )
  
  if( !isHPGe && lowres_shrink_roi_for_multipeaks && (initialfitpeaks.size() > 1 || lowstatregion) )
    initialfitpeaks = lowres_shrink_roi( initialfitpeaks, dataH, lowstatregion, automated );
  
  if( isHPGe && initialfitpeaks.size()==1 )
    initialfitpeaks = highres_shrink_roi( initialfitpeaks, dataH, lowstatregion, automated );
  
  if( initialfitpeaks.size() == 1 )
  {
    std::shared_ptr<const PeakDef> p = initialfitpeaks[0];
    
    bool passed;
    if( isHPGe )
      passed = check_highres_single_peak_fit( p, dataH, automated );
    else
      passed = check_lowres_single_peak_fit( p, dataH, lowstatregion, automated );
    
    if( passed )
      return std::make_pair( initialfitpeaks, coFitPeaks );
    
    return pair<PeakShrdVec,PeakShrdVec>();
  }//if( !isHPGe && fitpeaks.size()==1 )
  
  
  pair<PeakShrdVec,PeakShrdVec> answer;
  answer.first = initialfitpeaks;
  answer.second = coFitPeaks;
  
  return answer;
}//searchForPeakFromUser(...)


void secondDerivativePeakCanidates( const std::shared_ptr<const Measurement> data,
                                   const bool isHPGe,
                                   size_t start_channel,
                                   size_t end_channel,
                                   std::vector< std::tuple<float,float,float> > &results )
{
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  //If we're debuging things, lets make sure we printout information from just
  //  one function call contiguously to the log.
  static std::mutex s_fcnmutex;
  std::lock_guard<std::mutex> lock( s_fcnmutex );
  
  ofstream ouputfile( "secondDerivativePeakCanidates.txt", ios::app );
  //  ostream &debugstrm = cout;
  ostream &debugstrm = ouputfile;
#endif
  
  results.clear();
  
  if( !data->num_gamma_channels() )
    return;
  
  const size_t nchannel = data->num_gamma_channels();
  
  if( start_channel >= nchannel )
    start_channel = 0;
  if( end_channel <= start_channel || end_channel >= (nchannel-1) )
    end_channel = nchannel - 2;
  
  const double threshold_FOM = 1.3;
  //const double range_nsigma_thresh = isHPGe ? 5.0 : 2.75;
  const float pos_sum_threshold_sf = -0.01f;
  //We will let one bin fluctate negative to avoid fluctations near threshold_FOM
  //Untested for HPGe data.
  //  Currently (20141209) for low res data, this should be kept the same as
  //  in find_roi_for_2nd_deriv_candidate(...) {although all this is under
  //  development}.
  const size_t nFluxuate = isHPGe ? 2 : 2;
  assert( nFluxuate >= 1 );
  
  
  //XXX: using middle energy range so peak finding will always be consistent
  //     despite intput range
  //     Should keep consistent with find_roi_for_2nd_deriv_candidate()
  const size_t midbin = data->num_gamma_channels() / 2;// (start_channel + end_channel) / 2;
  const float midenergy = data->gamma_channel_center( midbin );
  const float midbinwidth = data->gamma_channel_width( midbin );
  
  const int order = isHPGe ? 3 : 2;
  const size_t side_bins = isHPGe ? 4 : std::max( size_t(5), static_cast<size_t>( 0.022f*midenergy/midbinwidth + 0.5f ) );
  
  
  vector<float> second_deriv;
  smoothSpectrum( data, static_cast<int>(side_bins), order, 2, second_deriv );
  
  //Since the peak resolution changes a lot for low-resolution spectra, if
  //  side_bins is to large, it will wipe out low energy features, meaning we
  //  wont detect low energy peaks
  //Note this code should be kept the same as in
  //  find_roi_for_2nd_deriv_candidate(...) while all of this is in development
  if( !isHPGe && side_bins > 5 && nchannel >= 512
     && start_channel < (nchannel/15) )
  {
    //We should also have a minbimum statistics requirment here.
    const size_t index = std::max( (nchannel/15), side_bins );
    vector<float> second_deriv_lower;
    smoothSpectrum( data, 4, order, 2, second_deriv_lower );
    
    for( size_t i = 0; i < (index-side_bins); ++i )
      second_deriv[i] = second_deriv_lower[i];
    
    //transition over 'side_bins' between the smoothings.
    for( size_t i = 0; i < side_bins; ++i )
    {
      const float factor = float(i+1) / float(side_bins+1);
      const size_t current = index - side_bins + i;
      second_deriv[current] = factor*second_deriv[current]
      + (1.0f-factor)*second_deriv_lower[current];
    }
    
  }//if( !isHPGe )
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  {
    static std::mutex secondderivfilelock;
    std::lock_guard<std::mutex> lock( secondderivfilelock );
    ofstream secondderiv( "secondderiv.csv" );
    secondderiv << "Energy,Counts" << endl;
    const size_t nchannel = data->gamma_counts()->size();
    for( size_t bin = 0; bin < nchannel; ++bin )
      secondderiv << data->gamma_channel_lower(bin) << "," << second_deriv[bin] << endl;
    vector<float> smoothed;
    smoothSpectrum( data, side_bins, order, 0, smoothed );
    ofstream smoothedfile( "smoothed.csv" );
    smoothedfile << "Energy,Counts" << endl;
    for( size_t bin = 0; bin < nchannel; ++bin )
      smoothedfile << data->gamma_channel_lower(bin) << "," << smoothed[bin] << endl;
    DebugLog(debugstrm) << "Made smoothed.csv and secondderiv.csv" << "\n";
    //    cerr << "side_bins=" << side_bins << ", order=" << order << endl;
  }
#endif
  
  //XXX: the below 1.5 is empiracally found, I'm not entirely sure where
  //     comes from...and infact might be higher
  const double amp_fake_factor = 1.5;
  
  
  const vector<float> &energies = *data->gamma_channel_energies();
  
  size_t minbin = 0, firstzero = 0, secondzero = 0;
  float secondsum = 0.0f, minval = 9999999999.9f;
  
  for( size_t channel = start_channel; channel <= end_channel; ++channel )
  {
    const float secondDeriv = second_deriv[channel]; //Not dividing by binwidth^2 here,
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 3 )
    DebugLog(debugstrm) << data->gamma_channel_lower(channel) << ": secondDeriv=" << secondDeriv
    << ", secondsum=" << secondsum << ", minval=" << minval << ", firstzero="
    << firstzero << ", channel=" << channel << "\n";
#endif
    
    bool secondSumPositive = true;
    float positivesum = 0.0f;
    for( size_t i = 0; i < nFluxuate; ++i )
    {
      if( (channel+i) <= end_channel )
      {
        const bool above = (second_deriv[channel+i] > 0.0f);
        if( above )
          positivesum += second_deriv[channel+i];
        secondSumPositive &= above;
      }
    }
    
    //Rather than using 'pos_sum_threshold_sf*secondsum', it might be better to
    //  use something invlolving sqrt(secondsum) since this makes a bit more
    //  sense.
    //Also, positivesum can probably also thresholded off of some sort of more
    //  absolute quantity
    secondSumPositive &= (positivesum > pos_sum_threshold_sf*secondsum);
    
    if( secondSumPositive && (minval < 99999999999.9)
       && (secondsum!=0.0) && (firstzero>0)
       && ((channel-firstzero)>2) )
    {
      secondzero = channel;
      
      const double mean = data->gamma_channel_center(minbin);
      const double sigma = 0.5*(data->gamma_channel_center(secondzero)
                                - data->gamma_channel_center(firstzero));
      
      const double deriv_sigma = 0.5*( secondzero - firstzero );
      const double part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                       boost::math::constants::e<double>() ) )
      / ( deriv_sigma * deriv_sigma );
      const double amplitude = -amp_fake_factor * secondsum / part;
      
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 2 )
      DebugLog(debugstrm) << "firstzero=" << data->gamma_channel_center(firstzero)
      << ", secondzero=" << data->gamma_channel_center(secondzero) << "\n";
#endif
      
      const float lowerEnengy = static_cast<float>( mean - 3.0*sigma );
      const float upperEnergy = static_cast<float>( mean + 3.0*sigma );
      
      double data_area = data->gamma_integral( lowerEnengy, upperEnergy );
      double est_sigma = sqrt( std::max(data_area,1.0) );
      
      //In principle we would want to use the (true) continuums area to derive
      //  the est_sigma from, but for practical purposes I think this can give
      //  us false positives fairly often
      
      const double figure_of_merit = 0.68*amplitude / est_sigma;
      
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 2 )
      //double cont_area = peak->offset_integral( lowerEnengy, upperEnergy );
      //double new_amp = data_area - cont_area;
      //DebugLog(debugstrm) << "mean=" << mean << ", amplitude=" << amplitude
      //<< ", sigma=" << peak->sigma()
      //<< ", lowerEnengy=" << lowerEnengy << ", upperEnergy=" << upperEnergy
      //<< ", cont_area=" << cont_area << ", data_area=" << data_area
      //<< ", new_amp=" << new_amp << ", FOM=" << figure_of_merit << "\n"
      //<< "\n";
#endif
      
      if( figure_of_merit > threshold_FOM )
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 2 )
        DebugLog(debugstrm) << "Candidate a peak with mean=" << mean
        << ", width=" << sigma
        << ", amplitude=" << amplitude
        << ", new_amp=" << amplitude
        << ", ROI start=" << lowerEnengy << ", ROI end=" << upperEnergy
        << ", figure_of_merit=" << figure_of_merit
        << "\n\tKeeping" << "\n";
#endif
        
        bool passescuts = true;
        if( !isHPGe && energies[minbin] < 130.0f )
        {
          //look 2 sigma forward and make sure the data has dropped enough
          const size_t p2sigmbin = data->find_gamma_channel( mean + 1.5*sigma );
          const float p2binlower = data->gamma_channel_lower( p2sigmbin );
          const float p2binupper = data->gamma_channel_upper( p2sigmbin );
          
          const double p2gausheight = amplitude*PeakDists::gaussian_integral( mean, sigma, p2binlower, p2binupper );
          const float p2contents = data->gamma_channel_content( p2sigmbin );
          
          const float meanbinlower = data->gamma_channel_lower( minbin );
          const float meanbinupper = data->gamma_channel_upper( minbin );
          const double meangausheight = amplitude*PeakDists::gaussian_integral( mean, sigma, meanbinlower, meanbinupper );
          const float meancontents = data->gamma_channel_content( minbin );
          const double expecteddiff = meangausheight - p2gausheight;
          const float actualdiff = meancontents - p2contents;
          
          //We dont expect the continuum to be changing rate of any more than
          //  1/2 the rate of peaks (this rate made up from my head, so take
          //  it with a grain of salt).  Should incorporate stat uncertainites
          //  here as well!
          passescuts = (actualdiff > 0.45*expecteddiff);
          
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          if( !passescuts )
            DebugLog(debugstrm) << "Failing candidate peak, mean: " << mean
            << ", meancontents=" << meancontents
            << ", p2contents=" << p2contents << ", actualdiff="
            << actualdiff << ", expecteddiff=" << expecteddiff << "\n";
#endif
        }//if( !isHPGe && energies[minbin] < 130.0f )
        
        if( !isHPGe && passescuts )
        {
          //For really wide peaks, we'll let there be a bit more fluxuation,
          //  since sometimes
          size_t nflux = max( nFluxuate, ((secondzero-firstzero)/side_bins) + 1 );
          
          //Look forward and backwards to see if the sum of the second
          //  derivative, while positive, is roughly comparable to the negative
          //  sum.  We expect the positive-summs on either side of the negative
          //  region to add up to about the same area as the negative region.
          float nextpositivesum = 0.0;
          for( size_t i = channel; i <= end_channel; ++i )
          {
            bool secondSumNegative = true;
            for( size_t j = 0; j < nflux && ((i+j) < end_channel); ++j )
              secondSumNegative &= (second_deriv[i+j] < 0.0f);
            if( secondSumNegative )
              break;
            nextpositivesum += second_deriv[i];
          }//for( size_t i = 0; i <= end_channel; ++i )
          
          
          float prevpositivesum = 0.0;
          for( size_t i = firstzero; i > 0; --i )
          {
            bool secondSumNegative = true;
            for( size_t j = 0; j < nflux && ((i-j) > 0); ++j )
              secondSumNegative &= (second_deriv[i-j] < 0.0f);
            if( secondSumNegative )
              break;
            prevpositivesum += second_deriv[i];
          }//for( size_t i = 0; i <= end_channel; ++i )
          
          //If the current candidate peak is a result of compton backscatter,
          //  then it is likely the next positive region will be large, while
          //  the previous positive region will be small (we expect them to be
          //  about equal, and sum to be equal in magnitide to the negative
          //  region).
          const float nextratio = (-nextpositivesum/secondsum);
          const float prevratio = (-prevpositivesum/secondsum);
          passescuts = (nextratio < 4.0 || prevratio > 0.2)
          && ((nextratio+prevratio)>0.3);
          
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          if( !passescuts )
            DebugLog(debugstrm) << "Failing candidate peak, mean=" << mean
            << ", prevsumratios=" << (-prevpositivesum/secondsum)
            << ", possumratios=" << (-nextpositivesum/secondsum)
            << ", ((secondzero-firstzero)/side_bins)=" << ((secondzero-firstzero)/side_bins)
            << "\n";
#endif
        }//if( !isHPGe )
        
        if( passescuts )
          results.push_back( std::tuple<float,float,float>{mean, sigma, amplitude} );
      }//if( region we were just in passed_threshold )
      
      secondsum = 0.0;
      minval = 9999999999.9f;
      minbin = secondzero = firstzero = 0;
    }else
    {
      bool belowzero = true, goingnegative = true, abovezero = true;
      for( size_t i = 0; i < nFluxuate; ++i )
      {
        if( (channel+i+1) < nchannel )
          goingnegative &= (second_deriv[channel+i+1] < 0.0f);
        if( channel >= i )
        {
          belowzero &= (second_deriv[channel-i] <= 0.0f);
          abovezero &= (second_deriv[channel-i] > 0.0f);
        }
      }//for( size_t i = 0; i < nFluxuate; ++i )
      
      if( channel /*&& (firstzero==0)*/ && !firstzero && goingnegative )
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
}//secondDerivativePeakCanidates(...)


std::vector<std::shared_ptr<PeakDef> > secondDerivativePeakCanidatesWithROI( std::shared_ptr<const Measurement> dataH,
                                                                            const bool isHPGe,
                                                          size_t start_channel,
                                                          size_t end_channel )
{
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  //If we're debuging things, lets make sure we printout information from just
  //  one function call contiguously to the log.
  static std::mutex s_fcnmutex;
  std::lock_guard<std::mutex> lock( s_fcnmutex );
  
  ofstream ouputfile( "secondDerivativePeakCanidatesWithROI.txt", ios::app );
//  ostream &debugstrm = cout;
  ostream &debugstrm = ouputfile;
#endif
  
  vector<std::shared_ptr<PeakDef> > candidates;
  
  if( !dataH->num_gamma_channels() )
    return candidates;
  
  const size_t nchannel = dataH->num_gamma_channels();
  
  if( start_channel >= nchannel )
    start_channel = 0;
  if( end_channel <= start_channel || end_channel >= (nchannel-1) )
    end_channel = nchannel - 2;
  
  const double threshold_FOM = 1.3;
  const double range_nsigma_thresh = isHPGe ? 5.0 : 2.75;
  const float pos_sum_threshold_sf = -0.01f;
  //We will let one bin fluctate negative to avoid fluctations near threshold_FOM
  //Untested for HPGe data.
  //  Currently (20141209) for low res data, this should be kept the same as
  //  in find_roi_for_2nd_deriv_candidate(...) {although all this is under
  //  development}.
  const size_t nFluxuate = isHPGe ? 2 : 2;
  assert( nFluxuate >= 1 );
  
  
  //XXX: using middle energy range so peak finding will always be consistent
  //     despite intput range
  //     Should keep consistent with find_roi_for_2nd_deriv_candidate()
  const size_t midbin = dataH->num_gamma_channels() / 2;// (start_channel + end_channel) / 2;
  const float midenergy = dataH->gamma_channel_center( midbin );
  const float midbinwidth = dataH->gamma_channel_width( midbin );
  
  const int order = isHPGe ? 3 : 2;
  const size_t side_bins = isHPGe ? 4 : std::max( size_t(5), static_cast<size_t>( 0.022f*midenergy/midbinwidth + 0.5f ) );
  
  
  vector<float> second_deriv;
  smoothSpectrum( dataH, static_cast<int>(side_bins), order, 2, second_deriv );
  
  //Since the peak resolution changes a lot for low-resolution spectra, if
  //  side_bins is to large, it will wipe out low energy features, meaning we
  //  wont detect low energy peaks
  //Note this code should be kept the same as in
  //  find_roi_for_2nd_deriv_candidate(...) while all of this is in development
  if( !isHPGe && side_bins > 5 && nchannel >= 512
      && start_channel < (nchannel/15) )
  {
    //We should also have a minbimum statistics requirment here.
    const size_t index = std::max( (nchannel/15), side_bins );
    vector<float> second_deriv_lower;
    smoothSpectrum( dataH, 4, order, 2, second_deriv_lower );
    
    for( size_t i = 0; i < (index-side_bins); ++i )
      second_deriv[i] = second_deriv_lower[i];
    
    //transition over 'side_bins' between the smoothings.
    for( size_t i = 0; i < side_bins; ++i )
    {
      const float factor = float(i+1) / float(side_bins+1);
      const size_t current = index - side_bins + i;
      second_deriv[current] = factor*second_deriv[current]
      + (1.0f-factor)*second_deriv_lower[current];
    }
    
  }//if( !isHPGe )
  
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
  {
    static std::mutex secondderivfilelock;
    std::lock_guard<std::mutex> lock( secondderivfilelock );
    ofstream secondderiv( "secondderiv.csv" );
    secondderiv << "Energy,Counts" << endl;
    const size_t nchannel = dataH->gamma_counts()->size();
    for( size_t bin = 0; bin < nchannel; ++bin )
      secondderiv << dataH->gamma_channel_lower(bin) << "," << second_deriv[bin] << endl;
    vector<float> smoothed;
    smoothSpectrum( dataH, side_bins, order, 0, smoothed );
    ofstream smoothedfile( "smoothed.csv" );
    smoothedfile << "Energy,Counts" << endl;
    for( size_t bin = 0; bin < nchannel; ++bin )
      smoothedfile << dataH->gamma_channel_lower(bin) << "," << smoothed[bin] << endl;
    DebugLog(debugstrm) << "Made smoothed.csv and secondderiv.csv" << "\n";
//    cerr << "side_bins=" << side_bins << ", order=" << order << endl;
  }
#endif
  
  //XXX: the below 1.5 is empiracally found, I'm not entirely sure where
  //     comes from...and infact might be higher
  const double amp_fake_factor = 1.5;
  

  const vector<float> &energies = *dataH->gamma_channel_energies();
  
  size_t minbin = 0, firstzero = 0, secondzero = 0;
  float secondsum = 0.0f, minval = 9999999999.9f;
  
  for( size_t channel = start_channel; channel <= end_channel; ++channel )
  {
    const float secondDeriv = second_deriv[channel]; //Not dividing by binwidth^2 here,
    
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 3 )
        DebugLog(debugstrm) << dataH->gamma_channel_lower(channel) << ": secondDeriv=" << secondDeriv
        << ", secondsum=" << secondsum << ", minval=" << minval << ", firstzero="
        << firstzero << ", channel=" << channel << "\n";
#endif
    
    bool secondSumPositive = true;
    float positivesum = 0.0f;
    for( size_t i = 0; i < nFluxuate; ++i )
    {
      if( (channel+i) <= end_channel )
      {
        const bool above = (second_deriv[channel+i] > 0.0f);
        if( above )
          positivesum += second_deriv[channel+i];
        secondSumPositive &= above;
      }
    }
    
    //Rather than using 'pos_sum_threshold_sf*secondsum', it might be better to
    //  use something invlolving sqrt(secondsum) since this makes a bit more
    //  sense.
    //Also, positivesum can probably also thresholded off of some sort of more
    //  absolute quantity
    secondSumPositive &= (positivesum > pos_sum_threshold_sf*secondsum);
    
    if( secondSumPositive && (minval < 99999999999.9)
       && (secondsum!=0.0) && (firstzero>0)
       && ((channel-firstzero)>2) )
    {
      secondzero = channel;
      
      const double mean = dataH->gamma_channel_center(minbin);
      const double sigma = 0.5*(dataH->gamma_channel_center(secondzero)
                                - dataH->gamma_channel_center(firstzero));
      
      const double deriv_sigma = 0.5*( secondzero - firstzero );
      const double part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                       boost::math::constants::e<double>() ) )
      / ( deriv_sigma * deriv_sigma );
      const double amplitude = -amp_fake_factor * secondsum / part;
      
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 2 )
      DebugLog(debugstrm) << "firstzero=" << dataH->gamma_channel_center(firstzero)
      << ", secondzero=" << dataH->gamma_channel_center(secondzero) << "\n";
#endif
      
      std::shared_ptr<PeakDef> peak( new PeakDef( mean, sigma, amplitude ) );
      
      double lowerEnengy, upperEnergy;
      findROIEnergyLimits( lowerEnengy, upperEnergy, *peak, dataH, isHPGe );
      
      //Clamp the ROI to not get ridiculous
      lowerEnengy = std::max( lowerEnengy, mean-5.0*sigma );
      upperEnergy = std::min( upperEnergy, mean+5.0*sigma );
      
      std::shared_ptr<PeakContinuum> continuum = peak->continuum();
      continuum->calc_linear_continuum_eqn( dataH, mean, lowerEnengy, upperEnergy, 2, 2 );
      
      const size_t lowchannel = dataH->find_gamma_channel( lowerEnengy );
      const size_t highchannel = dataH->find_gamma_channel( upperEnergy );
      lowerEnengy = dataH->gamma_channel_lower( lowchannel );
      upperEnergy = dataH->gamma_channel_upper( highchannel );
      
      double data_area = dataH->gamma_channels_sum( lowchannel, highchannel );
      double est_sigma = sqrt( std::max(data_area,1.0) );
      
      //In principle we would want to use the (true) continuums area to derive
      //  the est_sigma from, but for practical purposes I think this can give
      //  us false positives fairly often
      
      const double figure_of_merit = 0.68*peak->amplitude()/est_sigma;

#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 2 )
      double cont_area = peak->offset_integral( lowerEnengy, upperEnergy, dataH );
      double new_amp = data_area - cont_area;
      DebugLog(debugstrm) << "mean=" << mean << ", amplitude=" << amplitude
            << ", sigma=" << peak->sigma()
            << ", lowerEnengy=" << lowerEnengy << ", upperEnergy=" << upperEnergy
            << ", cont_area=" << cont_area << ", data_area=" << data_area
            << ", new_amp=" << new_amp << ", FOM=" << figure_of_merit << "\n"
            << "\n";
#endif
      
      assert( fabs(sigma-peak->sigma()) < FLT_EPSILON );
      bool rangeOk = ((upperEnergy-lowerEnengy)>range_nsigma_thresh*sigma);
      
      //A problem here is if we are near another definite peak, then the ROI
      //  for this one may be small, causing the check for this below to fail.
      if( !rangeOk && isHPGe && candidates.size() )
      {
        const PeakDef &last = *candidates.back();
        const double s = 0.5*(last.sigma() + sigma);
        const double dx = mean - last.mean();
        rangeOk = ( ((dx/s) < range_nsigma_thresh) && ((dx/s) > 1.75));
      }//if( !rangeOk && candidates.size() )
      
      
      if( (figure_of_merit > threshold_FOM) && rangeOk )
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 2 )
        DebugLog(debugstrm) << "Candidate a peak with mean=" << peak->mean()
                << ", width=" << peak->sigma()
                << ", amplitude=" << peak->amplitude()
                << ", new_amp=" << new_amp
                << ", ROI start=" << lowerEnengy << ", ROI end=" << upperEnergy
                << ", cont_area=" << cont_area
                << ", figure_of_merit=" << figure_of_merit
                 << "\n\tKeeping" << "\n";
#endif
        
        bool passescuts = true;
        if( !isHPGe && energies[minbin] < 130.0f )
        {
          //look 2 sigma forward and make sure the data has dropped enough
          const size_t p2sigmbin = dataH->find_gamma_channel( mean + 1.5*sigma );
          const float p2binlower = dataH->gamma_channel_lower( p2sigmbin );
          const float p2binupper = dataH->gamma_channel_upper( p2sigmbin );
          const double p2gausheight = peak->gauss_integral( p2binlower, p2binupper );
          const float p2contents = dataH->gamma_channel_content( p2sigmbin );
          
          const float meanbinlower = dataH->gamma_channel_lower( minbin );
          const float meanbinupper = dataH->gamma_channel_upper( minbin );
          const double meangausheight = peak->gauss_integral( meanbinlower, meanbinupper );
          const float meancontents = dataH->gamma_channel_content( minbin );
          const double expecteddiff = meangausheight - p2gausheight;
          const float actualdiff = meancontents - p2contents;
          
          //We dont expect the continuum to be changing rate of any more than
          //  1/2 the rate of peaks (this rate made up from my head, so take
          //  it with a grain of salt).  Should incorporate stat uncertainites
          //  here as well!
          passescuts = (actualdiff > 0.45*expecteddiff);
          
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          if( !passescuts )
            DebugLog(debugstrm) << "Failing candidate peak, mean: " << mean
                 << ", meancontents=" << meancontents
                 << ", p2contents=" << p2contents << ", actualdiff="
                 << actualdiff << ", expecteddiff=" << expecteddiff << "\n";
#endif
        }//if( !isHPGe && energies[minbin] < 130.0f )
        
        if( !isHPGe && passescuts )
        {
          //For really wide peaks, we'll let there be a bit more fluxuation,
          //  since sometimes 
          size_t nflux = max( nFluxuate, ((secondzero-firstzero)/side_bins) + 1 );
          
          //Look forward and backwards to see if the sum of the second
          //  derivative, while positive, is roughly comparable to the negative
          //  sum.  We expect the positive-summs on either side of the negative
          //  region to add up to about the same area as the negative region.
          float nextpositivesum = 0.0;
          for( size_t i = channel; i <= end_channel; ++i )
          {
            bool secondSumNegative = true;
            for( size_t j = 0; j < nflux && ((i+j) < end_channel); ++j )
              secondSumNegative &= (second_deriv[i+j] < 0.0f);
            if( secondSumNegative )
              break;
            nextpositivesum += second_deriv[i];
          }//for( size_t i = 0; i <= end_channel; ++i )
          
          
          float prevpositivesum = 0.0;
          for( size_t i = firstzero; i > 0; --i )
          {
            bool secondSumNegative = true;
            for( size_t j = 0; j < nflux && ((i-j) > 0); ++j )
              secondSumNegative &= (second_deriv[i-j] < 0.0f);
            if( secondSumNegative )
              break;
            prevpositivesum += second_deriv[i];
          }//for( size_t i = 0; i <= end_channel; ++i )
          
          //If the current candidate peak is a result of compton backscatter,
          //  then it is likely the next positive region will be large, while
          //  the previous positive region will be small (we expect them to be
          //  about equal, and sum to be equal in magnitide to the negative
          //  region).
          const float nextratio = (-nextpositivesum/secondsum);
          const float prevratio = (-prevpositivesum/secondsum);
          passescuts = (nextratio < 4.0 || prevratio > 0.2)
                       && ((nextratio+prevratio)>0.3);

#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
          if( !passescuts )
            DebugLog(debugstrm) << "Failing candidate peak, mean=" << peak->mean()
                 << ", prevsumratios=" << (-prevpositivesum/secondsum)
                 << ", possumratios=" << (-nextpositivesum/secondsum)
                 << ", ((secondzero-firstzero)/side_bins)=" << ((secondzero-firstzero)/side_bins)
                 << "\n";
#endif
        }//if( !isHPGe )
        
        if( passescuts )
          candidates.push_back( peak );
      }//if( region we were just in passed_threshold )
      
      secondsum = 0.0;
      minval = 9999999999.9f;
      minbin = secondzero = firstzero = 0;
    }else
    {
      bool belowzero = true, goingnegative = true, abovezero = true;
      for( size_t i = 0; i < nFluxuate; ++i )
      {
        if( (channel+i+1) < nchannel )
          goingnegative &= (second_deriv[channel+i+1] < 0.0f);
        if( channel >= i )
        {
          belowzero &= (second_deriv[channel-i] <= 0.0f);
          abovezero &= (second_deriv[channel-i] > 0.0f);
        }
      }//for( size_t i = 0; i < nFluxuate; ++i )
      
      if( channel /*&& (firstzero==0)*/ && !firstzero && goingnegative )
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
  
  return candidates;
}//secondDerivativePeakCanidatesWithROI( std::shared_ptr<const Measurement> dataH,



void get_chi2_and_dof_for_roi( double &chi2, double &dof,
                              const std::shared_ptr<const SpecUtils::Measurement> &data,
                              const vector<PeakDef *> &peaks )
{
  dof = chi2 = 0.0;
  if( peaks.empty() || !data || !data->channel_energies() || data->channel_energies()->empty() )
    return;
  
  assert( peaks[0] );
  const std::shared_ptr<const PeakContinuum> continuum = peaks[0]->continuum();
  
  for( size_t i = 1; i < peaks.size(); ++i )
  {
    assert( peaks[i] );
    assert( continuum == peaks[i]->continuum() );
  }
  
  const double lx = continuum->lowerEnergy();
  const double ux = continuum->upperEnergy();
  
  const size_t startchannel = data->find_gamma_channel( lx + 0.0000001 );
  const size_t endchannel = data->find_gamma_channel( ux - 0.0000001 );
  const size_t numchannel = (endchannel >= startchannel) ? (1 + endchannel - startchannel) : size_t(0);
  const vector<float> &energies = *data->channel_energies();
  assert( startchannel < energies.size() );
  assert( endchannel < energies.size() );
  
  vector<double> gauss_counts( std::max(numchannel,size_t(0)), 0.0 );
  for( size_t i = 0; i < peaks.size(); ++i )
    peaks[i]->gauss_integral( &(energies[startchannel]), &(gauss_counts[0]), numchannel );
  
  for( size_t i = 0; i < numchannel; ++i )
  {
    const size_t channel = startchannel + i;
    const double xbinlow = energies[channel];
    const double xbinup  = energies[channel+1];
    double nfitpeak = gauss_counts[i];
    
    const double ndata = data->gamma_channel_content( channel );
    const double ncontinuim = continuum->offset_integral( xbinlow, xbinup, data );
    
    const double datauncert = std::max( ndata, 1.0 );
    const double nabove = (ndata - ncontinuim - nfitpeak);
    chi2 += nabove*nabove / datauncert;
  }//for( size_t i = 0; i < numchannel; ++i )
  

  int nfitsigma = 0, nfitamp = 0, nfitmean = 0;
  for( auto p : peaks )
  {
    nfitsigma += p->fitFor(PeakDef::CoefficientType::Sigma);
    nfitmean += p->fitFor(PeakDef::CoefficientType::Mean);
    nfitamp += p->fitFor(PeakDef::CoefficientType::GaussAmplitude);
  }

  dof = 1.0*(endchannel - startchannel);// - 2.0 - 1.0*(peaks.size()-1);// - 1.0*continuum->type();
  
  //We fit the amplitudes independantly
  dof -= nfitamp;

  //The means are kinda, sorta, independent
  dof -= nfitmean;
  
  //We fit with the sigmas tied together, with a max of two parameters to control across the whole
  //  ROI - actually its a bit more complicated than this, but good enough for right now
  dof -= std::min( nfitsigma, 2 );
  
  for( const bool fit : continuum->fitForParameter() )
    dof -= (fit ? 1.0 : 0.0);
  
  dof -= 1.0;
  
  if( dof < 1.0 )
    dof = 1.0;
  
  /// \TODO: (20200718) more accurately calculate the DOF
  
  //Note: Previous to 20200718, we used the following formula for DOF... not quite sure why:
  //  1.0*(endchannel - startchannel) - 2.0 - 1.0*(peaks.size()-1); - 1.0*continuum->type();
}//get_chi2_and_dof_for_roi(...)



double set_chi2_dof( std::shared_ptr<const Measurement> data,
                    std::vector<PeakDef> &fitpeaks,
                    const size_t startpeak, const size_t npeaks )
{
  double totalDOF = 0;
  //It would be nice to use chi2_for_region(...) to actually compute the chi2
  //  for the region
  
  map< std::shared_ptr<const PeakContinuum>, vector<PeakDef *> > roigroups;
  for( size_t i = startpeak; i < npeaks && i < fitpeaks.size(); ++i )
  {
    if( fitpeaks[i].gausPeak() )
      roigroups[fitpeaks[i].continuum()].push_back( &(fitpeaks[i]) );
  }
  
  for( auto i = begin(roigroups); i != end(roigroups); ++i )
  {
    const vector<PeakDef *> &peakptrs = i->second;
    assert( peakptrs.size() );
    
    double chi2, dof;
    get_chi2_and_dof_for_roi( chi2, dof, data, peakptrs );
    
    totalDOF += dof;
    const double chi2Dof = chi2 / dof;
    
    for( PeakDef * const peak : peakptrs )
    {
      peak->set_coefficient( chi2Dof, PeakDef::Chi2DOF );
      peak->set_uncertainty( 0.0, PeakDef::Chi2DOF );
    }
    
    //cout << "Set chi2dof=" << chi2Dof << endl;
  }//for( loop over ROIs )
  
  return totalDOF;
}//void set_chi2_dof( )


void fitPeaks( const std::vector<PeakDef> &all_near_peaks,
              const double stat_threshold,
              const double hypothesis_threshold,
              std::shared_ptr<const Measurement> data,
              std::vector<PeakDef> &fitpeaks,
              const std::vector<PeakDef> &all_fixedpeaks,
              bool amplitudeOnly,
              const bool isHPGe ) throw()
{
  try
  {
    fitpeaks.clear();
      
    //We have to separate out non-gaussian peaks since they cant enter the
    //  fitting methods
    vector<PeakDef> fixedpeaks, near_peaks, datadefined_peaks;
    
    std::shared_ptr<const Measurement> continuum;
    
    fixedpeaks.reserve( all_fixedpeaks.size() );
    near_peaks.reserve( all_near_peaks.size() );
    
    for( const PeakDef &p : all_near_peaks )
    {
      if( p.gausPeak() )
        near_peaks.push_back( p );
      else
        datadefined_peaks.push_back( p );
      
      if( !!p.continuum()->externalContinuum() )
        continuum = p.continuum()->externalContinuum();
    }//for( const PeakDef &p : all_near_peaks )
    
    for( const PeakDef &p : all_fixedpeaks )
    {
      if( p.gausPeak() )
        fixedpeaks.push_back( p );
      else
        datadefined_peaks.push_back( p );
      
      if( !!p.continuum()->externalContinuum() )
        continuum = p.continuum()->externalContinuum();
    }//for( const PeakDef &p : all_fixedpeaks )
    
    if( near_peaks.empty() )
      return;
    
    unique_copy_continuum( near_peaks );
    unique_copy_continuum( fixedpeaks );  //prob not necassary, but JIC
      
    //Need to make sure near_peaks and fixedpeaks are all gaussian (if not
    //  seperate them out, and add them in later).  If fitpeaks is non-gaussian
    //  ignore it or throw an exception.
    
    double lowx( 0.0 ), highx( 0.0 );
    
    {
      const PeakDef &lowgaus = near_peaks.front();
      const PeakDef &highgaus = near_peaks.back();
      
      double dummy = 0.0;
      findROIEnergyLimits( lowx, dummy, lowgaus, data, isHPGe );
      findROIEnergyLimits( dummy, highx, highgaus, data, isHPGe );
    }
    
    const int npeaks = static_cast<int>( near_peaks.size() + fixedpeaks.size() );
    PeakFitChi2Fcn chi2Fcn( npeaks, data, continuum );
    chi2Fcn.useReducedChi2( false );
    
    //    ROOT::Minuit2::Minuit2Minimizer fitter( ROOT::Minuit2::kMigrad );
    //    fitter.SetFunction( chi2Fcn );
    //    fitter.SetMaxFunctionCalls( 2500 );
    //    fitter.SetTolerance( 0.1 );
    //    fitter.SetPrintLevel( -1 );
    
    /*
     static std::mutex smutex;
     {
     std::lock_guard<std::mutex> lock( smutex );
     cerr << "fitPeaks(): we have " << near_peaks.size() << " near peaks, and "
     << fixedpeaks.size() << " fixedpeaks\nNearPeaks:" << endl;
     for( const PeakDef &peak : near_peaks )
       cerr << "\t" << peak.mean() << endl;
     cerr << "FixedPeaks:" << endl;
     for( const PeakDef &peak : fixedpeaks )
       cerr << "\t" << peak.mean() << endl;
     }
     */
    
    ROOT::Minuit2::MnUserParameters inputPrams;
    PeakFitChi2Fcn::AddPeaksToFitterMethod method = (amplitudeOnly
                                                     ? PeakFitChi2Fcn::kRefitPeakParameters
                                                     : PeakFitChi2Fcn::kFitForPeakParameters);
    
    PeakFitChi2Fcn::addPeaksToFitter( inputPrams, near_peaks, data, method, isHPGe );
    PeakFitChi2Fcn::addPeaksToFitter( inputPrams, fixedpeaks, data,
                                     PeakFitChi2Fcn::kFixPeakParameters, isHPGe );
    
    if( inputPrams.VariableParameters() == 0 )
    {
      fitpeaks = near_peaks;
      return;
    }//if( inputPrams.VariableParameters() == 0 )
    
    
    ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
    ROOT::Minuit2::MnStrategy strategy( 2 ); //0 low, 1 medium, >=2 high
    ROOT::Minuit2::MnMinimize fitter( chi2Fcn, inputParamState, strategy );
    //fitter.SetPrecision(<#double#>)
    
    unsigned int maxFcnCall = 5000;
    double tolerance = 2.5;
    tolerance = 0.5;
    ROOT::Minuit2::FunctionMinimum minimum = fitter( maxFcnCall, tolerance );
    const ROOT::Minuit2::MnUserParameters &fitParams = fitter.Parameters();
    //  minimum.IsValid()
    //      ROOT::Minuit2::MinimumState minState = minimum.State();
    //      ROOT::Minuit2::MinimumParameters minParams = minState.Parameters();
    
    //    cerr << endl << endl << "EDM=" << minimum.Edm() << endl;
    //    cerr << "MinValue=" <<  minimum.Fval() << endl << endl;
    
    if( !minimum.IsValid() )
      minimum = fitter( maxFcnCall, tolerance );
      
    if( !minimum.IsValid() )
    {
      //XXX - should we try to re-fit here? Or do something to handle the
      //      failure in some reasonable way?
#ifndef NDEBUG
      cerr << endl << endl << "status is not valid"
           << "\n\tHasMadePosDefCovar: " << minimum.HasMadePosDefCovar()
           << "\n\tHasAccurateCovar: " << minimum.HasAccurateCovar()
           << "\n\tHasReachedCallLimit: " << minimum.HasReachedCallLimit()
           << "\n\tHasValidCovariance: " << minimum.HasValidCovariance()
           << "\n\tHasValidParameters: " << minimum.HasValidParameters()
           << "\n\tIsAboveMaxEdm: " << minimum.IsAboveMaxEdm()
           << endl;
      if( minimum.IsAboveMaxEdm() )
        cout << "\t\tEDM=" << minimum.Edm() << endl;
#endif
    }//if( !minimum.IsValid() )
    
    
    vector<double> fitpars = fitParams.Params();
    vector<double> fiterrors = fitParams.Errors();
    chi2Fcn.parametersToPeaks( fitpeaks, &fitpars[0], &fiterrors[0] );
      
    //Lets try to keep whether or not to fit parameters should be the same for
    //  the output peaks as the input peaks.
    //Note that this doesnt account for peaks swapping with eachother in the fit
    if( fitpeaks.size() != (near_peaks.size()+fixedpeaks.size()) )
      throw std::logic_error( "fitPeaks(...): unexpected number of peaks" );
        
    for( size_t i = 0; i < near_peaks.size(); ++i )
      fitpeaks[i].inheritUserSelectedOptions( near_peaks[i], true );
    for( size_t i = 0; i < fixedpeaks.size(); ++i )
      fitpeaks[i+near_peaks.size()].inheritUserSelectedOptions( fixedpeaks[i], true );
      
    set_chi2_dof( data, fitpeaks, 0, near_peaks.size() );
      
    std::sort( fitpeaks.begin(), fitpeaks.end(), &PeakDef::lessThanByMean );
        
    vector<PeakDef> *all_peaks = &fitpeaks;
    std::unique_ptr< vector<PeakDef> > all_peaks_ptr;
    if( fixedpeaks.size() )
    {
      all_peaks = new vector<PeakDef>( fitpeaks );
      all_peaks_ptr.reset( all_peaks );
      all_peaks->insert( all_peaks->end(), fixedpeaks.begin(), fixedpeaks.end() );
      std::sort( all_peaks->begin(), all_peaks->end(), &PeakDef::lessThanByMean );
    }//if( fixedpeaks.size() )
    
    for( size_t i = 1; i <= fitpeaks.size(); ++i ) //Note weird convntion of index
    {
      const PeakDef *peak = &(fitpeaks[i-1]);
      
      // Dont remove peaks whos amplitudes we arent fitting
      if( !peak->fitFor(PeakDef::GaussAmplitude) )
        continue;
      
      // Dont enforce a significance test if we are refitting the peak, and
      //  Sigma and Mean are fixed - the user probably knows what they
      //  are are doing.
      if( (method == PeakFitChi2Fcn::kRefitPeakParameters)
        && !peak->fitFor(PeakDef::Sigma)
        && !peak->fitFor(PeakDef::Mean) )
      {
        continue;
      }

      const bool significant = chi2_significance_test( *peak, stat_threshold, hypothesis_threshold,
                                                          *all_peaks, data );
      if( !significant )
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cerr) << "\tPeak at mean=" << peak->mean()
                       << "is being discarded for not being significant"
                       << "\n";
#endif
        fitpeaks.erase( fitpeaks.begin() + --i );
      }//if( !significant )
    }//for( size_t i = 1; i < fitpeaks.size(); ++i )
        
    bool removed_peak = false;
    for( size_t i = 1; i < fitpeaks.size(); ++i ) //Note weird convention of index
    {
      PeakDef *this_peak = &(fitpeaks[i-1]);
      PeakDef *next_peak = &(fitpeaks[i+1-1]);
          
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
          fitpeaks.erase( fitpeaks.begin() + i - 1 );
        else
          fitpeaks.erase( fitpeaks.begin() + i );
            
        i = i - 1; //incase we have multiple close peaks in a row
      }//if( (mean_diff/min_sigma) < 1.0 ) / else
    }//for( size_t i = 1; i < fitpeaks.size(); ++i )
        
    if( removed_peak )
      fitPeaks( fitpeaks, stat_threshold, hypothesis_threshold,
                data, fitpeaks, fixedpeaks, false, isHPGe );
        
    if( datadefined_peaks.size() )
    {
      fitpeaks.insert( fitpeaks.end(),
                      datadefined_peaks.begin(), datadefined_peaks.end() );
      std::sort( fitpeaks.begin(), fitpeaks.end(), &PeakDef::lessThanByMean );
    }
        
    return;
  }catch( std::exception &e )
  {
    cerr << "fitPeaks(...)\n\tSerious programming logic error: caught"
         << " exception where I really shouldnt have.  what()=" << e.what()
         << endl;
  }catch(...)
  {
    cerr << "fitPeaks(...)\n\tSerious programming logic error: caught"
         << " unknown exception where I really shouldnt have." << endl;
  }//try/catch()
        
  //We will only reach here if there was no exception, so since never expect
  //  this to actually happen, just assign the results to be same as the input
  fitpeaks = all_near_peaks;
}//vector<PeakDef> fitPeaks(...);
        
        
        
std::vector<PeakDef> peaksInRange( const double lowx,
                                   const double highx,
                                   const double nsigma,
                                   const std::vector<PeakDef> &inputs )
{
  vector<PeakDef> answer;
          
  for( const PeakDef &peak : inputs )
  {
    const double peakmin  = (peak.gausPeak() ? (peak.mean() - nsigma*peak.sigma()) : peak.lowerX());
    const double peakmax  = (peak.gausPeak() ? (peak.mean() + nsigma*peak.sigma()) : peak.upperX());
            
            //1   ----------        x
            //          ++++++++    peak
            //
            //2       ----------    x
            //    ++++++++          peak
            //
            //3   ----------------  x
            //        ++++++        peak
            //
            //4         -------     x
            //      +++++++++++++++ peak
            
    if( (peakmin<=highx) && (peakmax>=lowx) )  //covers case 1 and 3
      answer.push_back( peak );
    else if( (peakmin<=lowx) && (peakmax>=lowx) ) //covers case 2
      answer.push_back( peak );
  }//for( const PeakDef &peak : inputs )
          
  std::sort( answer.begin(), answer.end(), &PeakDef::lessThanByMean );
          
  return answer;
}//peaksInRange(...)
        
        
        std::vector<PeakDef> peaksTouchingRange( double lowx, double highx,
                                                const std::vector<PeakDef> &inputs )
        {
          if( highx < lowx )
            std::swap( lowx, highx );
          
          vector<PeakDef> answer;
          
          for( const PeakDef &peak : inputs )
          {
            const double peakmin  = peak.lowerX();
            const double peakmax  = peak.upperX();
            
            //1   ----------        x
            //          ++++++++    peak
            //
            //2       ----------    x
            //    ++++++++          peak
            //
            //3   ----------------  x
            //        ++++++        peak
            //
            //4         -------     x
            //      +++++++++++++++ peak
            
            if( (peakmin<=highx) && (peakmax>=lowx) )  //covers case 1 and 3
              answer.push_back( peak );
            else if( (peakmin<=lowx) && (peakmax>=lowx) ) //covers case 2
              answer.push_back( peak );
          }//for( const PeakDef &peak : inputs )
          
          std::sort( answer.begin(), answer.end(), &PeakDef::lessThanByMean );
          
          return answer;
        }//peaksTouchingRange(...)
        
        
PeakShrdVec peaksTouchingRange( double lowx, double highx,
                                const PeakShrdVec &inputs )
{
  typedef std::shared_ptr<const PeakDef> PeakDefShrdPtr;
  
  if( highx < lowx )
    std::swap( lowx, highx );
          
  PeakShrdVec answer;
          
  for( const PeakDefShrdPtr &peak : inputs )
  {
    const double peakmin  = peak->lowerX();
    const double peakmax  = peak->upperX();
            
    if( (peakmin<=lowx) && (peakmax>=highx) )
      answer.push_back( peak );
    else if( (peak->mean()>=lowx) && (peak->mean()<=highx) )
      answer.push_back( peak );
    else if( (peakmin>=lowx) && (peakmin<=highx) )
      answer.push_back( peak );
    else if( (peakmax>=lowx) && (peakmax<=highx) )
      answer.push_back( peak );
  }//for( const PeakDef &peak : inputs )
          
  std::sort( answer.begin(), answer.end(), &PeakDef::lessThanByMeanShrdPtr );
          
  return answer;
}//peaksTouchingRange(...)
        
        
double evaluate_polynomial( const double x,
                                 const std::vector<double> &poly_coeffs )
{
  double y = 0.0;
  for( size_t i = 0; i < poly_coeffs.size(); ++i )
    y += poly_coeffs[i] * std::pow( x, double(i) );
  return y;
}//
      
      
double fit_to_polynomial( const float *x, const float *data, const size_t nbin,
                                 const int polynomial_order,
                                 std::vector<double> &poly_coeffs,
                                 std::vector<double> &coeff_uncerts )
{
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //Implementation is quite inneficient
  using namespace boost::numeric;
  const int poly_terms = polynomial_order + 1;
  ublas::matrix<double> A( nbin, poly_terms );
  ublas::vector<double> b( nbin );
  
  for( size_t row = 0; row < nbin; ++row )
  {
    const double uncert = (data[row] > MIN_CHANNEL_UNCERT ? sqrt( data[row] ) : 1.0);
    b(row) = (data[row] > 0.0 ? sqrt( data[row] ) : 0.0);
    for( int col = 0; col < poly_terms; ++col )
      A(row,col) = std::pow( double(x[row]), double(col)) / uncert;
  }//for( int col = 0; col < poly_terms; ++col )
  
  const ublas::matrix<double> A_transpose = ublas::trans( A );
  const ublas::matrix<double> alpha = prod( A_transpose, A );
  ublas::matrix<double> C( alpha.size1(), alpha.size2() );
  const bool success = matrix_invert( alpha, C );
  if( !success )
    throw runtime_error( "fit_to_polynomial(...): trouble inverting matrix" );
  
  const ublas::vector<double> beta = prod( A_transpose, b );
  const ublas::vector<double> a = prod( C, beta );
  
  poly_coeffs.resize( poly_terms );
  coeff_uncerts.resize( poly_terms );
  for( int coef = 0; coef < poly_terms; ++coef )
  {
    poly_coeffs[coef] = a(coef);
    coeff_uncerts[coef] = std::sqrt( C(coef,coef) );
  }//for( int coef = 0; coef < poly_terms; ++coef )
  
  double chi2 = 0;
  for( size_t bin = 0; bin < nbin; ++bin )
  {
    double y_pred = 0.0;
    for( int i = 0; i < poly_terms; ++i )
      y_pred += a(i) * std::pow( double(x[bin]), double(i) );
    const double uncert = (data[bin] > MIN_CHANNEL_UNCERT ? sqrt( data[bin] ) : 1.0);
    chi2 += std::pow( (y_pred - data[bin]) / uncert, 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2;
}//double fit_to_polynomial(...)
        
        
double fit_amp_and_offset( const float *x, const float *data, const size_t nbin,
                                  const int num_polynomial_terms,
                                  const bool step_continuum,
                                  const double ref_energy,
                                  const vector<double> &means,
                                  const vector<double> &sigmas,
                                  const vector<PeakDef> &fixedAmpPeaks,
                                  const PeakDef::SkewType skew_type,
                                  const double *skew_parameters,
                                  std::vector<double> &amplitudes,
                                  std::vector<double> &continuum_coeffs,
                                  std::vector<double> &amplitudes_uncerts,
                                  std::vector<double> &continuum_coeffs_uncerts )
{
  // TODO: Need to switch to using Eigen::SVD for this function - it is much more stable and predictable
  //       See RelActCalcAuto::fit_continuum(...) for example of using this.
  if( sigmas.size() != means.size() )
    throw runtime_error( "fit_amp_and_offset: invalid input" );
  
  if( step_continuum && ((num_polynomial_terms < 2) || (num_polynomial_terms > 4)) )
    throw runtime_error( "fit_amp_and_offset: Only 2 to 4 terms are supported for step continuums" );
  
  if( num_polynomial_terms < 0 )
    throw runtime_error( "fit_amp_and_offset: continuum must have at least 0 (e.g., no continuum) terms" );
  
  if( num_polynomial_terms > 4 )
    throw runtime_error( "fit_amp_and_offset: you asked for a higher order polynomial continuum than reasonable" );
  
  assert( (skew_type == PeakDef::SkewType::NoSkew) || skew_parameters );
  if( !skew_parameters && (skew_type != PeakDef::SkewType::NoSkew) )
    throw std::logic_error( "Skew pars not provided" );
  
  //Using variable names of section 15.4 of Numerical Recipes, 3rd edition
  //
  //Implementation is quite inefficient
  //
  //Current implementation is not necessarily numerically the most accurate or
  //  the best when there are near degeneracies.  Should switch to using SVD for
  //  solving.
  //
  //When uncertainties are not needed (e.g., while fitting means), could save
  //  maybe a factor of three or so in computations.
  //
  //Could use Eigen and do things a bit better, see
  // https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  //
  
  
//For experimentation to see how much all the allocations are actually slowing
//  things done - horribly thread un-safe, and will explode if used from GUI!
#define TRIAL_MINIMIZE_MATRIX_ALLOCATIONS 0
  
  
  using namespace boost::numeric;
  const size_t npeaks = sigmas.size();
  const size_t npoly = static_cast<size_t>( num_polynomial_terms );
  const int nfit_terms = static_cast<int>( npoly + npeaks );
  
#if( TRIAL_MINIMIZE_MATRIX_ALLOCATIONS )
#warning "TRIAL_MINIMIZE_MATRIX_ALLOCATIONS Is enabled - this will blow up if making an InterSpec GUI"
  static ublas::matrix<double> A( nbin, nfit_terms );
  static ublas::vector<double> b( nbin );
  
  if( A.size1() != nbin || A.size2() != nfit_terms )
    A.resize( nbin, nfit_terms, false );
  if( b.size() != nbin )
    b.resize( nbin, false );
#else
  ublas::matrix<double> A( nbin, nfit_terms );
  ublas::vector<double> b( nbin );
#endif
  
  //  cerr << endl << "Input: " << ref_energy << ", ";
  //  for( size_t i = 0; i < npeaks; ++i )
  //    cerr << "{" << means[i] << ", " << sigmas[i] << "}, ";
  //  cerr << endl << endl;
  
  double roi_data_sum = 0.0, step_cumulative_data = 0.0;
  
  const double roi_lower = x[0];
  const double roi_upper = x[nbin];
  
  for( size_t row = 0; row < nbin; ++row )
    roi_data_sum += std::max( data[row], 0.0f );
  
  const double avrg_data_val = roi_data_sum / nbin;
  
  // For the RelACtAuto calcs, there can easily be 100 peaks, so for this case we'll calculate
  //  their contribution multithreaded.
  //  I have no idea if this is the optimal way to do the calculation, but better than not using
  //  threads at all.
  //  If we dont have many fixed peaks, we'll avoid the extra allocation (totally unchecked if
  //  this actually saves anything perceptible)
  //  We'll divide the channels into nthread ranges, and then inside each thread,
  //  compute its range for each peak - this way we can cut the number of calls to the `erf`
  //  function in half by calling a more optimized version of the peak integral function
  //
  const size_t nfixedpeak = fixedAmpPeaks.size();
  const bool do_mt_fixed_peak = (nfixedpeak > 0); // always use this method if fixed peaks
  vector<double> mt_fixed_peak_contrib( do_mt_fixed_peak ? nbin : size_t(0), 0.0f );
  
  
  if( do_mt_fixed_peak )
  {
    // 20250127: it looks like calling `pool.join()` is causing significant and unreasonable delays
    //           on macOS (using GCD, at least).  This is likely a problem with
    //           `SpecUtilsAsync::ThreadPool` - I would guess when creating ThreadPools inside of
    //           other ThreadPools, but for the moment will just do this single threaded, which is
    //           like 20 times faster for an example problem
    
    //SpecUtilsAsync::ThreadPool pool;
    double * const fixed_contrib = &(mt_fixed_peak_contrib[0]);
    
    /*
    // 20240325: for complicated problems, it kinda looks like multiple cores arent being used that
    //           efficiently (perhaps because each thread is only getting <10 channels to work on).
    //           Perhaps it would be better to put each peak into a thread, for all channels, then
    //           sum results at the end (using vector operations).
    const size_t nthread = std::thread::hardware_concurrency();
    const size_t nbin_per_thread = 1 + (nbin / nthread);
    
    for( size_t start_channel = 0; start_channel < nbin; start_channel += nbin_per_thread )
    {
      const float *this_x = x + start_channel;
      double *this_contribs = fixed_contrib + start_channel;
      const size_t end_channel = start_channel + nbin_per_thread;
      const size_t this_nchannel = std::min( end_channel, nbin ) - start_channel; //make sure we dont go out of range
      
      assert( (this_nchannel > 0) && (this_nchannel <= nbin_per_thread) );
      
      pool.post( [this_x, this_contribs, this_nchannel, &fixedAmpPeaks](){
        for( size_t peak_index = 0; peak_index < fixedAmpPeaks.size(); ++peak_index )
        {
          const PeakDef &peak = fixedAmpPeaks[peak_index];
          peak.gauss_integral( this_x, this_contribs, this_nchannel );
        }
      } );
    }//for( size_t start_channel = 0; start_channel < nbin; start_channel += nbin_per_thread )
    */
    
    //std::mutex result_mutex;
    for( size_t peak_index = 0; peak_index < fixedAmpPeaks.size(); ++peak_index )
    {
      /*
      pool.post( [&result_mutex, peak_index, &fixedAmpPeaks, nbin, x, fixed_contrib](){
        vector<double> this_contribs(nbin, 0.0);
        fixedAmpPeaks[peak_index].gauss_integral( x, &(this_contribs[0]), nbin );
        
        std::lock_guard<std::mutex> lock( result_mutex );
        for( size_t i = 0; i < nbin; ++i )
          fixed_contrib[i] += this_contribs[i];
      });
       */
      
      fixedAmpPeaks[peak_index].gauss_integral( x, &(fixed_contrib[0]), nbin );
    }//for( size_t peak_index = 0; peak_index < fixedAmpPeaks.size(); ++peak_index )
    
    
    //pool.join();
  }//if( do_mt_fixed_peak )
  
  
  for( size_t row = 0; row < nbin; ++row )
  {
    double dataval = data[row];
    
    const double x0 = x[row];
    const double x1 = x[row+1];
    
    const double x0_rel = x0 - ref_energy;
    const double x1_rel = x1 - ref_energy;
    
    //const double uncert = (dataval > 0.0 ? sqrt(dataval) : 1.0);
    // If we are background subtracting a spectrum, we can end up with bins with really
    //  small values, like 0.0007, which, even one of would mess the whole fit up if
    //  we take its uncertainty to be its square-root, so in this case we will, fairly arbitrarily
    //  we want to use an uncert of 1.
    //  However, there are also highly scaled spectra, whose all values are really small, so
    //  in this case we want to do something more reasonable.
    // TODO: evaluate these choices of thresholds and tradeoffs, more better
    double uncert = (dataval > MIN_CHANNEL_UNCERT ? sqrt(dataval) : 1.0);
  
    if( step_continuum )
      step_cumulative_data += dataval;
    
/*
    // If data is zero, or negative, lets look for the nearest non-zero bin, within 5 bins of here
    //  This situation might happen more often after a hard background subtraction.
    //
    //  TODO: this doesnt fix the one-example I was looking at - perhaps need to try ignoring a
    //        region. Maybe try looking for regions that are anomalously low (e.g., essentially
    //        zero), and just ignore them if they are extremely below surrounding region.
    //
    const double non_poisson_threshold = 0.5; //FLT_EPSILON
    const double significant_stats_threshold = 5.0;
  
    if( dataval < non_poisson_threshold )
    {
      double nearest_data = dataval;
      
      for( size_t i = 1; (i <= 5) && (nearest_data < significant_stats_threshold); ++i )
      {
        if( ((row + i) < nbin) && (data[row+i] > significant_stats_threshold) )
          nearest_data = data[row+i];
        
        if( i <= row && (data[row-i] > significant_stats_threshold) )
          nearest_data = data[row-i];
      }//for( size_t i = 1; (i <= 5) && (nearest_data < non_poisson_threshold); ++i )
      
      if( nearest_data > non_poisson_threshold )
        uncert = sqrt(nearest_data);
    }//if( dataval < FLT_EPSILON )
*/
    
    if( do_mt_fixed_peak )
    {
      assert( mt_fixed_peak_contrib.size() == nbin );
      dataval -= mt_fixed_peak_contrib[row];
    }else if( !fixedAmpPeaks.empty() )
    {
      for( size_t i = 0; i < fixedAmpPeaks.size(); ++i )
      {
        // See multithreaded implementation above for reasoning; the logic and width should
        //  match above.
        const PeakDef &peak = fixedAmpPeaks[i];
        const double mean = peak.mean();
        const double integration_width = 8*peak.sigma();
        if( (x1 >= (mean - integration_width)) && (x0 <= (mean + integration_width)) )
          dataval -= fixedAmpPeaks[i].gauss_integral( x0, x1 );
      }
    }//if( do_mt_fixed_peak )
    
    b(row) = ((dataval > 0.0 ? dataval : 0.0) / uncert);
    
    for( size_t col = 0; col < npoly; ++col )
    {
      const double exp = col + 1.0;
      
      if( step_continuum
          && ((num_polynomial_terms == 2) || (num_polynomial_terms == 3))
          && (col == (num_polynomial_terms - 1)) )
      {
        // This logic mirrors that of PeakContinuum::offset_integral(...), and code
        // If you change it in one place - change it in here, below, and in offset_integral.
        const double frac_data = (step_cumulative_data - 0.5*data[row]) / roi_data_sum;
        const double contribution = frac_data * (x1 - x0);
        
        A(row,col) = contribution / uncert;
      }else if( step_continuum && (num_polynomial_terms == 4) )
      {
        const double frac_data = (step_cumulative_data - 0.5*data[row]) / roi_data_sum;

        double contrib = 0.0;
        switch( col )
        {
          case 0: contrib = (1.0 - frac_data) * (x1_rel - x0_rel);                     break;
          case 1: contrib = 0.5 * (1.0 - frac_data) * (x1_rel*x1_rel - x0_rel*x0_rel); break;
          case 2: contrib = frac_data * (x1_rel - x0_rel);                             break;
          case 3: contrib = 0.5 * frac_data * (x1_rel*x1_rel - x0_rel*x0_rel);         break;
          default: assert( 0 ); break;
        }//switch( col )
        
        A(row,col) = contrib / uncert;
      }else
      {
        const double contribution = (1.0/exp) * (pow(x1_rel,exp) - pow(x0_rel,exp));
        
        A(row,col) = contribution / uncert;
      }
    }//for( int order = 0; order < maxorder; ++order )
  }//for( size_t row = 0; row < nbin; ++row )
  
  
  // If we have more than 2 peaks (arbitrarily chosen), we'll compute the peak contributions
  //  in parallel.
  //  TODO: investigate performance impact of computing peak integrals multithread - e.g., should we do this only if a skew is being used?  Or if we have more than X peaks, etc.
  vector<double> unit_peak_counts( nbin * npeaks, 0.0 );
  
  const bool parallelize_peak_sum = (npeaks > 2);
  if( parallelize_peak_sum )
  {
    SpecUtilsAsync::ThreadPool pool;
    for( size_t i = 0; i < npeaks; ++i )
    {
      double *peak_areas = &(unit_peak_counts[i*nbin]);
      const double mean = means[i];
      const double sigma = sigmas[i];
      pool.post( [peak_areas,i,mean,sigma,skew_type,skew_parameters,nbin,x](){
        PeakDists::photopeak_function_integral( mean, sigma, 1.0, skew_type, skew_parameters,
                                             nbin, x, peak_areas );
      } );
    }//for( size_t i = 0; i < npeaks; ++i )
    pool.join();
  }else
  {
    for( size_t i = 0; i < npeaks; ++i )
    {
      double *peak_areas = &(unit_peak_counts[i*nbin]);
      PeakDists::photopeak_function_integral( means[i], sigmas[i], 1.0,
                                           skew_type, skew_parameters,
                                           nbin, x, peak_areas );
    }
  }//if( npeaks > 2 ) / else
  
  for( size_t i = 0; i < npeaks; ++i )
  {
    double *peak_areas = &(unit_peak_counts[i*nbin]);
    for( size_t channel = 0; channel < nbin; ++channel )
    {
      const double dataval = data[channel];
      const double uncert = (dataval > MIN_CHANNEL_UNCERT ? sqrt(dataval) : 1.0);
      A(channel,npoly + i) = peak_areas[channel] / uncert;
    }//for( size_t channel = 0; channel < nbin; ++channel )
  }//for( size_t i = 0; i < npeaks; ++i )
  
#if( TRIAL_MINIMIZE_MATRIX_ALLOCATIONS )
  const ublas::matrix<double> A_transpose = ublas::trans( A );
  static ublas::matrix<double> alpha;
  
  if( alpha.size1() != A.size2() )
    alpha.resize( A.size2(), A.size2(), false );
  
  prod( A_transpose, A, alpha );
  
  static ublas::matrix<double> C;
  if( C.size1() != alpha.size1() || C.size2() != alpha.size2() )
    C.resize( alpha.size1(), alpha.size2(), false );
#else
  const ublas::matrix<double> A_transpose = ublas::trans( A );
  const ublas::matrix<double> alpha = prod( A_transpose, A );
  ublas::matrix<double> C( alpha.size1(), alpha.size2() );
#endif
  
  bool success = false;
  
  try
  {
    //See http://viennacl.sourceforge.net/doc/least-squares_8cpp-example.html#a8 for example
    //  of better solving things following the below commented out lines
    //typedef boost::numeric::ublas::matrix<double>              MatrixType;
    //typedef boost::numeric::ublas::vector<double>              VectorType;
    //boost::numeric::ublas::range ublas_range(0, 3);
    //boost::numeric::ublas::matrix_range<MatrixType> ublas_R(ublas_A, ublas_range, ublas_range);
    //boost::numeric::ublas::vector_range<VectorType> ublas_b2(ublas_b, ublas_range);
    //boost::numeric::ublas::inplace_solve(ublas_R, ublas_b2, boost::numeric::ublas::upper_tag());

#if( TRIAL_MINIMIZE_MATRIX_ALLOCATIONS )
    static boost::numeric::ublas::matrix<double> local_alpha;
    if( local_alpha.size1() != alpha.size1() || local_alpha.size2() != alpha.size2() )
      local_alpha.resize( alpha.size1(), alpha.size2());
    local_alpha.assign( alpha );
    
    ublas::permutation_matrix<std::size_t> pm( local_alpha.size1() );
    //if( pm.size() != local_alpha.size1() )
    //  pm.resize( local_alpha.size1(), false );
    
    const size_t res = boost::numeric::ublas::lu_factorize(local_alpha, pm);
    //const size_t res = boost::numeric::ublas::axpy_lu_factorize(local_alpha, pm);  //sometimes faster?
    success = (res == 0);
    if( success )
    {
      C.assign( boost::numeric::ublas::identity_matrix<double>( local_alpha.size1() ) );
      lu_substitute(local_alpha, pm, C);
    }
#else
    success = matrix_invert( alpha, C );
#endif
  }catch( std::exception &e )
  {
#ifndef NDEBUG
    cerr << "fit_amp_and_offset(...): caught: " << e.what() << endl;
    cerr << "For means = {";
    for( double m : means )
      cerr << m << ", ";
    cerr << "}, sigmas={";
    for( double m : sigmas )
      cerr << m << ", ";
    cerr << "}" << endl;
    
    printf( "b=" );
    for( size_t row = 0; row < b.size(); ++row )
      printf( "%.2f, ", b(row) );
    printf( "\n" );
    
    printf( "Alpha=\n" );
    for( size_t row = 0; row < alpha.size1(); ++row )
    {
      for( size_t col = 0; col < alpha.size2(); ++col )
        printf( "%12.2f, ", alpha(row,col) );
      printf( "\n" );
    }
    printf( "\n" );
    
    
    printf( "\nC=\n" );
    for( size_t row = 0; row < C.size1(); ++row )
    {
      for( size_t col = 0; col < C.size2(); ++col )
        printf( "%12.2f, ", C(row,col) );
      printf( "\n" );
    }
    printf( "\n\n\n" );
#endif //#ifndef NDEBUG
  }//try / catch
  
  if( !success )
  {
#ifndef NDEBUG
    cerr << "For means = {";
    for( double m : means )
      cerr << m << ", ";
    cerr << "}, sigmas={";
    for( double m : sigmas )
      cerr << m << ", ";
    cerr << "}" << endl;
#endif //#ifndef NDEBUG
    throw runtime_error( "fit_amp_and_offset(...): trouble inverting matrix" );
  }
  
  const ublas::vector<double> beta = prod( A_transpose, b );
  const ublas::vector<double> a = prod( C, beta );
  
  continuum_coeffs.resize( npoly );
  continuum_coeffs_uncerts.resize( npoly );
  for( size_t coef = 0; coef < npoly; ++coef )
  {
    continuum_coeffs[coef] = a(coef);
    continuum_coeffs_uncerts[coef] = std::sqrt( C(coef,coef) );
  }//for( int coef = 0; coef < poly_terms; ++coef )
  
  amplitudes.resize( npeaks );
  amplitudes_uncerts.resize( npeaks );
  
  for( size_t i = 0; i < npeaks; ++i )
  {
    const size_t coef = npoly + i;
    amplitudes[i] = a(coef);
    amplitudes_uncerts[i] = std::sqrt( C(coef,coef) );
  }//for( size_t i = 0; i < npeaks; ++i )
  
  double chi2 = 0;
  step_cumulative_data = 0.0;
  for( size_t bin = 0; bin < nbin; ++bin )
  {
    const double x0 = x[bin];
    const double x1 = x[bin+1];
    
    double dataval = data[bin];
    
    if( step_continuum )
      step_cumulative_data += dataval;
    
    //TODO: I havent actually reasoned through the algorithm to see if this is the
    //      correct way to subtract off fixed-amplitude peaks.
    for( size_t i = 0; i < fixedAmpPeaks.size(); ++i )
      dataval -= fixedAmpPeaks[i].gauss_integral( x0, x1 );
    
    double y_pred = 0.0;
    for( size_t col = 0; col < npoly; ++col )
    {
      const double exp = col + 1.0;
      const double x0_rel = x0 - ref_energy;
      const double x1_rel = x1 - ref_energy;
      
      if( step_continuum
         && ((num_polynomial_terms == 2) || (num_polynomial_terms == 3))
         && (col == (num_polynomial_terms - 1)) )
      {
        // This logic mirrors that of PeakContinuum::offset_integral(...) and above code in this
        //  function that defines the matrix, see above for comments
        const double frac_data = (step_cumulative_data - 0.5*data[bin]) / roi_data_sum;
        const double contribution = frac_data * (x1 - x0);
        
        y_pred += a(col)*contribution;
      }else if( step_continuum && (num_polynomial_terms == 4) )
      {
        // This logic mirrors that of PeakContinuum::offset_integral(...) and above code in this
        //  function that defines the matrix, see above for comments
        
        const double frac_data = (step_cumulative_data - 0.5*data[bin]) / roi_data_sum;
        
        double contrib = 0.0;
        switch( col )
        {
          case 0: contrib = (1.0 - frac_data) * (x1_rel - x0_rel);                     break;
          case 1: contrib = 0.5 * (1.0 - frac_data) * (x1_rel*x1_rel - x0_rel*x0_rel); break;
          case 2: contrib = frac_data * (x1_rel - x0_rel);                             break;
          case 3: contrib = 0.5 * frac_data * (x1_rel*x1_rel - x0_rel*x0_rel);         break;
          default: assert( 0 ); break;
        }//switch( col )
        
        y_pred += a(col) * contrib;
      }else
      {
        y_pred += a(col) * (1.0/exp) * (pow(x1_rel,exp) - pow(x0_rel,exp));
      }//if( step_continuum ) / else
    }//for( int order = 0; order < maxorder; ++order )
    
    if( y_pred < 0.0 )
      y_pred = 0.0;
    
    for( size_t i = 0; i < npeaks; ++i )
    {
      const size_t col = npoly + i;
      y_pred += a(col) * PeakDists::gaussian_integral( means[i], sigmas[i], x0, x1 );
    }
    
    for( size_t i = 0; i < fixedAmpPeaks.size(); ++i )
      y_pred += fixedAmpPeaks[i].gauss_integral( x0, x1 );
    
    //    cerr << "bin " << bin << " predicted " << y_pred << " data=" << data[bin] << endl;
    const double uncert = (data[bin] > MIN_CHANNEL_UNCERT ? sqrt( data[bin] ) : 1.0);
    chi2 += std::pow( (y_pred - data[bin]) / uncert, 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2;
}//double fit_amp_and_offset(...)
        
        
        
bool chi2_significance_test( PeakDef peak,
                             const double withoutPeakDSigma,
                             const double chi2ratioRequired,
                             std::vector<PeakDef> other_peaks,
                             std::shared_ptr<const Measurement> data )
{
  if( !peak.gausPeak() )
    throw runtime_error( "chi2_significance_test can only evaluate gaussian defined peaks" );
          
  //need to edit paramsForPeak and zeroPeakParams to evaluate the chi2 to a max
  //  of 2.5 sigma from the mean
  peak.makeUniqueNewContinuum();
  
  // We will limit the test range to +-2.5 FWHM of the peak mean to do the Chi2 test so we capture,
  //  only the relevant part of the peak.  But the step-continuum peaks we need to keep the original
  //  energy range, or else it will alter each channels continuum values (and we will fail peak fits
  //  we shouldnt).
  //
  //TODO 20210726: I dont really remember why I limited this check to +-2.5 FWHM from peak mean, but
  //  I'm guessing its to avoid failing long continuums where the tails are badly matched to the
  //  data, but user wants in order to force things.  Should investigate/document more thoroughly.
  bool isStepContinuum = false;
  switch( peak.continuum()->type() )
  {
    case PeakContinuum::NoOffset: case PeakContinuum::Constant:
    case PeakContinuum::Linear:   case PeakContinuum::Quadratic:
    case PeakContinuum::Cubic:    case PeakContinuum::External:
      isStepContinuum = false;
      break;
    
    case PeakContinuum::FlatStep:     case PeakContinuum::LinearStep:
    case PeakContinuum::BiLinearStep:
      isStepContinuum = true;
      break;
  }//switch( peak.continuum()->type() )
  
  
  if( !isStepContinuum )
  {
    double ux = peak.upperX();
    double lx = peak.lowerX();
    ux = std::min( ux, peak.mean() + 2.5*peak.sigma() );
    lx = std::max( lx, peak.mean() - 2.5*peak.sigma() );
    peak.continuum()->setRange( lx, ux );
  }//if( !isStepContinuum )
  
  //For section of code below here, please see notes in searchForPeaks()
  std::vector<PeakDef>::iterator pos = other_peaks.end();
  for( size_t i = 0; i < other_peaks.size(); ++i )
  {
    const PeakDef &other = other_peaks[i];
    
    if( fabs(other.mean() - peak.mean()) < 0.0001
        && fabs(other.sigma() - peak.sigma()) < 0.0001 )
    {
      pos = other_peaks.begin() + i;
      break;
    }
  }//for( size_t i = 0; i < other_peaks.size(); ++i )
          
  //  pos = std::find( other_peaks.begin(), other_peaks.end(), peak );
  if( pos != other_peaks.end() )
    other_peaks.erase( pos );
          
  //Some basic quality checks
  if( IsNan(peak.mean()) || IsInf(peak.mean()) )
    return false;
  
  if( peak.gausPeak() )
  {
    if( IsNan(peak.sigma()) || IsInf(peak.sigma()) )
      return false;
    
    if( IsNan(peak.amplitude()) || IsInf(peak.amplitude()) )
      return false;
    
    if( peak.fitFor(PeakDef::CoefficientType::GaussAmplitude) && ((peak.amplitude() <= 0.0)) )
      return false;
    
    if( peak.fitFor(PeakDef::CoefficientType::Sigma) && ((peak.sigma() <= 0.0)) )
      return false;
  }//if( peak.gausPeak() )
          
  vector<double> paramsForOtherPeaks, zeroPeakParams, paramsForPeak;
          
  {//start code block to create PeakFitChi2Fcn parameters for other peaks
    const bool isHPGe = false; //doesnt matter, since peaks already have ROI range defined
    ROOT::Minuit2::MnUserParameters mnparams;
    PeakFitChi2Fcn::addPeaksToFitter( mnparams, other_peaks,
                                    data, PeakFitChi2Fcn::kFixPeakParameters, isHPGe );
    paramsForOtherPeaks = mnparams.Params();
  }//end code block to create PeakFitChi2Fcn parameters for other peaks
          
  {//start code block to create PeakFitChi2Fcn parameters for zero amp test peak
    PeakDef zeropeak = peak;
    zeropeak.setAmplitude( 0.0 );
    
    assert( zeropeak.continuum()->energyRangeDefined() );
    const bool isHPGe = false; //doesnt matter, since peaks already have ROI range defined
    ROOT::Minuit2::MnUserParameters mnparams;
    PeakFitChi2Fcn::addPeaksToFitter( mnparams, vector<PeakDef>(1,zeropeak),
                                      data, PeakFitChi2Fcn::kFixPeakParameters, isHPGe);
    zeroPeakParams = mnparams.Params();
            
    //I dont think this next little bit is necessary, but I left in to be sure
    //  to be compatible w/ legacy code pre 20131230
    if( !zeropeak.continuum()->energyRangeDefined() )
    {
      //double xmin, xmax;
      //findROIEnergyLimits( xmin, xmax, peak, data );
      
      assert( peak.continuum()->energyRangeDefined() );
      const double xmin = peak.lowerX();
      const double xmax = peak.upperX();
      
      zeroPeakParams[PeakFitChi2Fcn::RangeStartEnergy] = xmin;
      zeroPeakParams[PeakFitChi2Fcn::RangeEndEnergy]   = xmax;
    }//if( !zeropeak.continuum()->energyRangeDefined() )
  }//end code block to create PeakFitChi2Fcn parameters for zero amp test peak
          
  {//start code block to create PeakFitChi2Fcn parameters for test peak
    const bool isHPGe = false; //Doesnt matter, since peaks already have ROI range defined
    ROOT::Minuit2::MnUserParameters mnparams;
    PeakFitChi2Fcn::addPeaksToFitter( mnparams, vector<PeakDef>(1,peak),
                                      data, PeakFitChi2Fcn::kFixPeakParameters, isHPGe );
    paramsForPeak = mnparams.Params();
  }//end code block to create PeakFitChi2Fcn parameters for test peak
          
  const int npeaks = static_cast<int>( other_peaks.size() + 1 );
          
  std::shared_ptr<const Measurement> contnuum;
  if( !!peak.continuum()->externalContinuum() )
  {
    contnuum = peak.continuum()->externalContinuum();
  }else
  {
    for( size_t i = 0; !contnuum && (i<other_peaks.size()); ++i )
    {
      if( !!other_peaks[i].continuum()->externalContinuum() )
      {
        contnuum = other_peaks[i].continuum()->externalContinuum();
        break;
      }
    }
  }//
          
          
  PeakFitChi2Fcn chi2fcn( npeaks, data, contnuum );
  chi2fcn.useReducedChi2( true );
  
  vector<double> params = paramsForOtherPeaks;
  params.insert( params.end(), zeroPeakParams.begin(), zeroPeakParams.end() );
  
  const double withoutGausChi2 = chi2fcn( &(params[0]) );
  
  params = paramsForOtherPeaks;
  params.insert( params.end(), paramsForPeak.begin(), paramsForPeak.end() );
  
  const double withGausChi2 = chi2fcn( &(params[0]) );
  const double chi2Ratio = withoutGausChi2 / withGausChi2;
  
  const bool noDeltaRequired = ((withoutPeakDSigma <= 0.0) && (chi2ratioRequired <= 0.0));
  bool noRatioRequired = noDeltaRequired;
          
  if( withoutGausChi2 < 5 )
    noRatioRequired = true;
  
  
          
  //Dont require the ratio test to apply if peaks share a continuum
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  for( const PeakDef &p : other_peaks )
    noRatioRequired |= (continuum == p.continuum());

  const double deltaChi2 = withoutGausChi2 - withGausChi2;
          
/*
  static std::mutex s_mutex;
  {
    std::lock_guard<std::mutex> lock( s_mutex );
    cerr << "chi2_significance_test(mean=" << peak.mean() << "): withGausChi2="
         << withGausChi2 << ", withoutGausChi2=" << withoutGausChi2 << endl
         << "chi2Ratio=" << chi2Ratio << " (requires " << chi2ratioRequired
         << "), deltaChi2=" << deltaChi2
         << " (required " << withoutPeakDSigma << ")"
         << ", noRatioRequired=" << noRatioRequired
         << ", peak.mean()=" << peak.mean()
         << ", peak.sigma()=" << peak.sigma()
         << " peak.amplitude()=" << peak.amplitude()
         << endl << endl;
  }
 */
  
  
  return ((noRatioRequired || (chi2Ratio >= chi2ratioRequired))
          && (noDeltaRequired || (deltaChi2 > withoutPeakDSigma)));
}//bool chi2_significance_test( ... 0
        
        
namespace ExperimentalPeakSearch
{

AutoPeakSearchChi2Fcn::AutoPeakSearchChi2Fcn( std::shared_ptr<const Measurement> data,
                                             const std::vector<PeakDef > &fixed_peaks,
                                             const bool isHPGe )
: ROOT::Minuit2::FCNBase(),
  m_inited( false ),
  m_isHPGe( isHPGe )
{
  if( !data )
    throw runtime_error( "AutoPeakSearchChi2Fcn: invalid input for construction" );
  
  m_meas = data;
  m_x = data->channel_energies();
  m_y = data->gamma_counts();
  
  m_fixed_peaks = fixed_peaks;
  
  m_side_bins           = m_isHPGe ? 7    : 10;
  m_smooth_order        = m_isHPGe ? 3    : 2;
  m_second_deriv_thresh = m_isHPGe ? -0.02 : 0.05;
  m_stat_thresh         = m_isHPGe ? 1.0   : 1.3;
  m_width_thresh        = 0.0;//m_isHPGe ? 3.5  : 4.5;
  
  m_min_chi2_dof_thresh = 1.5;
  m_min_gross_counts_sig_thresh = 2.0;
  m_nsigma_near_group = 3.0;
}//AutoPeakSearchChi2Fcn


size_t AutoPeakSearchChi2Fcn::lower_spectrum_channel() const { return m_lower_channel; }
size_t AutoPeakSearchChi2Fcn::upper_spectrum_channel() const { return m_upper_channel; }
size_t AutoPeakSearchChi2Fcn::num_initial_candidate_peaks() const { return m_num_peakstotal; }

void AutoPeakSearchChi2Fcn::second_derivative( const vector<float> &input, vector<float> &results )
{
  results.clear();
  SavitzyGolayCoeffs sgcoeffs( m_side_bins, m_side_bins, m_smooth_order, 2 );
  sgcoeffs.smooth( &(input[0]), static_cast<int>(input.size()), results );
}//void smoothSpectrum(...)



std::vector<PeakDef> AutoPeakSearchChi2Fcn::candidate_peaks( const vector<float> &energies,
                                                            const vector<float> &channel_counts )
{
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
  static int nwrites = 0;
  if( nwrites++ < 3 )
    cerr << "Writing candidate peak info to files" << endl;
#endif
  
  vector<PeakDef> candidates;
  
  const int nchannel = static_cast<int>( channel_counts.size() );
  
  vector<float> second_deriv;
  second_derivative( channel_counts, second_deriv );
  
  {
    vector<float> results;
    SavitzyGolayCoeffs sgcoeffs( m_side_bins, m_side_bins, m_smooth_order, 0 );
    sgcoeffs.smooth( &(channel_counts[0]), static_cast<int>(channel_counts.size()), results );
    
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
    ofstream fileout( "smoothout.csv" );
    fileout << "Energy,Counts" << endl;
    for( size_t i = 0; i < second_deriv.size(); ++i )
    fileout << energies[i] << "," << results[i] << endl;
#endif
  }
  
  {
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
    ofstream fileout( "secondout.csv" );
    fileout << "Energy,Counts" << endl;
    for( size_t i = 0; i < second_deriv.size(); ++i )
    fileout << energies[i] << "," << second_deriv[i] << endl;
#endif
  }
  
  int minbin = -1, firstzero = -1, secondzero=-1;
  float secondsum = 0.0f, minval = 9999999999.9f;
  
  for( int bin = static_cast<int>(m_lower_channel); bin < static_cast<int>(m_upper_channel); ++bin )
  {
    const float secondDeriv = second_deriv[bin]; //Not dividing by binwidth^2 here,
    
    if( secondDeriv > m_second_deriv_thresh*secondsum
       && (minval < 99999999999.9)
       && (secondsum!=0.0) && (firstzero>0)
       && ((bin-firstzero)>2) )
    {
      secondzero = bin;
      
      const double mean =  0.5*(energies[minbin] + energies[minbin+1]);
      const double sigma = 0.5*(energies[secondzero] - energies[firstzero]);
      
      //XXX: the below 1.5 is empiracally found, I'm not entirely sure where
      //     comes from...and infact might be higher
      const double amp_fake_factor = 1.5;
      const double deriv_sigma = 0.5*( secondzero - firstzero );
      const double part = sqrt( 2.0 / ( boost::math::constants::pi<double>() *
                                       boost::math::constants::e<double>() ) )
      / ( deriv_sigma * deriv_sigma );
      const double amplitude = -amp_fake_factor * secondsum / part;
      
      
      PeakDef peak( mean, sigma, amplitude );
      
      size_t lower_channel, upper_channel;
      double lowerEnengy = -999.9, upperEnergy = -999.9;
      
      //Set the ROI width here according to the second derivative.
      if( m_isHPGe )
      {
        int i, j;
        
        i = firstzero-1;
        while( i > 0 && second_deriv[i] > 0.0 )
          --i;
        
        i = i-1;
        while( i > 0 && second_deriv[i] < 0.0 )
          --i;
        
        j = i-1;
        while( j > 0 && second_deriv[j] > 0.0 )
          --j;
        
        lower_channel = ((i + j) / 2) + ((i+j)%2);
        
        i = secondzero+1;
        while( i < nchannel && second_deriv[i] > 0.0 )
          ++i;
        
        i = i+1;
        while( i < nchannel && second_deriv[i] < 0.0 )
          ++i;
        
        upper_channel = i-1;
        lowerEnengy = 0.5*(energies[lower_channel] + energies[lower_channel+1]);
        upperEnergy = 0.5*(energies[upper_channel] + energies[upper_channel+1]);
      }else
      {
        /*
         lowerEnengy = 0.5*(energies[firstzero-1] + energies[firstzero]);
         for( lowerbin = firstzero-1; lowerbin > 0 && second_deriv[lowerbin] > 0.0f; --lowerbin )
         lowerEnengy = 0.5*(energies[lowerbin] + energies[lowerbin+1]);
         ++lowerbin;
         
         upperEnergy = 0.5*(energies[secondzero+1] + energies[secondzero+2]);
         for( upperbin = secondzero+1; upperbin < (nchannel-1) && second_deriv[upperbin] > 0.0f; ++upperbin )
         upperEnergy = 0.5*(energies[upperbin] + energies[upperbin+1]);
         */
        
        findROIEnergyLimits( lowerEnengy, upperEnergy, peak, m_meas, m_isHPGe );
        lower_channel = m_meas->find_gamma_channel( lowerEnengy );
        upper_channel = m_meas->find_gamma_channel( upperEnergy );
      }//if( m_isHPGe ) / else
      
      //Clamp the ROI to not get riducuolouse
      //      lowerEnengy = std::max( lowerEnengy, mean-5.0*sigma );
      //      upperEnergy = std::min( upperEnergy, mean+5.0*sigma );
      
      std::shared_ptr<PeakContinuum> continuum = peak.continuum();
      continuum->setRange( lowerEnengy, upperEnergy );
      //      continuum->calc_linear_continuum_eqn( dataH, lowerEnengy, upperEnergy, 1 );
      
      double data_area = 0.0;
      assert( upper_channel < channel_counts.size() );
      for( size_t i = lower_channel; i <= upper_channel; ++i )
      data_area += channel_counts[i];
      
      double est_sigma = sqrt( std::max(data_area,1.0) );
      const double figure_of_merit = 0.68*peak.amplitude()/est_sigma;
      
      
      //        cout << "Found lowerEnengy=" << lowerEnengy
      //        << ", upperEnergy=" << upperEnergy
      //        << " for mean=" << mean << ", sigma=" << sigma
      //        << ", FOM="<<figure_of_merit << ", est_sigma=" << est_sigma
      //        << ", peak.amplitude()=" << peak.amplitude()
      //        << endl;
      
      if( figure_of_merit > m_stat_thresh
         && ((upperEnergy-lowerEnengy)>m_width_thresh*peak.sigma()) )
      {
        //          cout << "Found lowerEnengy=" << lowerEnengy
        //               << ", upperEnergy=" << upperEnergy
        //               << " for mean=" << mean << ", sigma=" << sigma << endl;
        
        candidates.push_back( peak );
      }//if( region we were just in passed_threshold )
      
      secondsum = 0.0;
      minval = 9999999999.9f;
      minbin = secondzero = firstzero = -1;
    }else
    {
      if( bin && (firstzero <= 0) && second_deriv[bin]>0.0 && second_deriv[bin+1]<0.0 )
      {
        firstzero = bin;
      }
      
      if( firstzero > 0 )
      {
        secondsum += secondDeriv;
        
        if( secondDeriv < minval )
        {
          minbin = bin;
          minval = secondDeriv;
        }
      }//if( firstzero > 0 )
    }//if( we are out of region of interest) / else( in region of )
  }//for( loop over bins )
  
  return candidates;
}//smoothingCandidatePeaks( std::shared_ptr<const Measurement> dataH,



//init(): returns if fitting can proced;
bool AutoPeakSearchChi2Fcn::init()
{
  if( m_inited )
    throw runtime_error( "AutoPeakSearchChi2Fcn: you cant call init twice" );
  
  ExperimentalPeakSearch::find_spectroscopic_extent( m_meas, m_lower_channel, m_upper_channel );
  
  
  if( m_fixed_peaks.empty() )
  {
    m_candidates = candidate_peaks( *m_x, *m_y );
  }else
  {
    vector<float> y = *m_y;
    for( size_t i = 0; i < (y.size()-1); ++i )
    {
      const float &lowere = (*m_x)[i];
      const float &uppere = (*m_x)[i+1];
      const float mide = 0.5f * (lowere + uppere);
      
      for( size_t j = 0; j < m_fixed_peaks.size(); ++j )
      {
        if( mide >= m_fixed_peaks[j].lowerX() && mide <= m_fixed_peaks[j].upperX() )
        {
          const double peakarea = m_fixed_peaks[j].gauss_integral( lowere, uppere );
          y[i] = std::max( y[i]-peakarea, 0.0 );
        }
      }//for( size_t j = 0; j < m_fixed_peaks.size(); ++j )
    }//for( size_t i = 0; i < (y.size()-1); ++i )
    
    std::vector<PeakDef> candidates = candidate_peaks( *m_x, y );
    
    
    //Now go through and remove all the candates within 1 sigma of an input put
    for( size_t i = 0; i < candidates.size(); ++i )
    {
      size_t nearest;
      double neareste = 99999999.9;
      for( size_t j = 0; j < m_fixed_peaks.size(); ++j )
      {
        const double de = fabs( m_fixed_peaks[j].mean() - candidates[i].mean() );
        if( de < neareste )
        {
          nearest = j;
          neareste = de;
        }
      }
      
      if( m_fixed_peaks.size() )
      {
        const PeakDef &p = m_fixed_peaks[nearest];
        if( fabs(p.mean() - candidates[i].mean()) > 1.0*p.sigma() )
          m_candidates.push_back( candidates[i] );
      }else
      {
        m_candidates.push_back( candidates[i] );
      }
    }//for( size_t i = 0; i < candidates.size(); ++i )
  }//if( m_fixed_peaks.empty() ) / else
  
  if( m_candidates.empty() )
    return false;
  
#if( WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL > 2 )
  for( const PeakDef &p : m_candidates )
  {
    cout << "Initial candidate peak { mean: " << p.mean()
    << ", sigma: " << p.sigma()
    << ", amplitude: " << p.amplitude()
    << ", ROI_lower: " << p.continuum()->lowerEnergy()
    << ", ROI_upper: " << p.continuum()->upperEnergy()
    << " }" << endl;
  }//for( const PeakDef &p : m_candidates )
#endif
  
  const size_t npeaks = m_fixed_peaks.size() + m_candidates.size();
  
  if( npeaks == 1 )
    m_resolution_type = Polynomial0thOrder;
  //    else
  //      m_resolution_type = SqrtEnergy;
  else if( npeaks < 4 )
    m_resolution_type = Polynomial1stOrder;
  else if( npeaks < 5 || !m_isHPGe )
    m_resolution_type = Polynomial2ndOrder;
  else
    m_resolution_type = Polynomial3rdOrder;
  
  //now we need to group peaks together into common ROI
  vector<PeakDef> allpeaks;
  allpeaks.insert( allpeaks.end(), m_candidates.begin(), m_candidates.end() );
  allpeaks.insert( allpeaks.end(), m_fixed_peaks.begin(), m_fixed_peaks.end() );
  std::sort( allpeaks.begin(), allpeaks.end(), &PeakDef::lessThanByMean );
  
  m_num_peakstotal = allpeaks.size();
  
  vector<bool> isfixed;
  vector<PeakDef> currentgroup;
  currentgroup.push_back( allpeaks[0] );
  
  
  std::vector<PeakDef>::const_iterator lb, ub;
  lb = std::lower_bound( m_fixed_peaks.begin(), m_fixed_peaks.end(), allpeaks[0], &PeakDef::lessThanByMean );
  ub = std::upper_bound( m_fixed_peaks.begin(), m_fixed_peaks.end(), allpeaks[0], &PeakDef::lessThanByMean );
  isfixed.push_back( (lb != ub) );
  
  size_t numpeaksinserted = 1;
  
  for( size_t i = 1; i < allpeaks.size(); ++i )
  {
    PeakDef p = allpeaks[i];
    
    lb = std::lower_bound( m_fixed_peaks.begin(), m_fixed_peaks.end(), p, &PeakDef::lessThanByMean );
    ub = std::upper_bound( m_fixed_peaks.begin(), m_fixed_peaks.end(), p, &PeakDef::lessThanByMean );
    
    const PeakDef &prev = currentgroup.back();
    
    const bool overlapping_roi = (p.lowerX() < prev.upperX());
    const bool isnear = (0.5*(p.mean() - prev.mean()) / (p.sigma() + prev.sigma())) < m_nsigma_near_group;
    
    if( overlapping_roi && isnear )
    {
      std::shared_ptr<PeakContinuum> cont = currentgroup.back().continuum();
      cont->setRange( cont->lowerEnergy(), std::max(p.upperX(),cont->upperEnergy()) );
      p.setContinuum( cont );
    }else
    {
      //        cout << "Adding group of " << isfixed.size() << " to fit for" << endl;
      m_grouped_isfixed.push_back( isfixed );
      m_grouped_candidates.push_back( currentgroup );
      
      std::shared_ptr<PeakContinuum> cont = currentgroup.back().continuum();
      
      //XXX - deciding the order of the continuum is purely a guess right now
      PeakContinuum::OffsetType type = PeakContinuum::Linear;
      if( m_isHPGe )
      {
        if( currentgroup.size() > 2 )
          type = PeakContinuum::Quadratic;
        if( currentgroup.size() > 4 )
          type = PeakContinuum::Cubic;
      }else
      {
        if( currentgroup.size() > 2 )
          type = PeakContinuum::Cubic;
      }
      
      cont->setType( type );
      
      isfixed.clear();
      currentgroup.clear();
    }//if( p.lowerX() < currentgroup.back().upperX() ) / else
    
    ++numpeaksinserted;
    isfixed.push_back( (lb != ub) );
    currentgroup.push_back( p );
  }//for( size_t i = 0; i < allpeaks.size(); ++i )
  
  assert( numpeaksinserted == npeaks );
  
  std::shared_ptr<PeakContinuum> cont = currentgroup.back().continuum();
  
  PeakContinuum::OffsetType type = PeakContinuum::Linear;
  if( m_isHPGe )
  {
    if( currentgroup.size() > 2 )
      type = PeakContinuum::Quadratic;
    if( currentgroup.size() > 4 )
      type = PeakContinuum::Cubic;
  }else
  {
    if( currentgroup.size() == 2 )
      type = PeakContinuum::Quadratic;
    if( currentgroup.size() > 2 )
      type = PeakContinuum::Cubic;
  }
  
  cont->setType( type );
  
  m_grouped_isfixed.push_back( isfixed );
  m_grouped_candidates.push_back( currentgroup );
  
  const float *x_begin = &(*m_x)[0];
  const float *x_end = (&(*m_x)[0]) + m_x->size();
  const float *y_begin = &(*m_y)[0];
  
  m_nbins_used = 0;
  for( size_t i = 0; i < m_grouped_candidates.size(); ++i )
  {
    std::shared_ptr<const PeakContinuum> cont = m_grouped_candidates[i][0].continuum();
    const double lowere = cont->lowerEnergy();
    const double uppere = cont->upperEnergy();
    
    const float *x_start = std::upper_bound( x_begin, x_end, lowere );
    if( x_start > x_begin )
      --x_start;
    const float *x_finish = std::upper_bound( x_begin, x_end, uppere );
    if( x_finish < x_end )
      ++x_finish;
    
    const size_t startbin = x_start - x_begin;
    const size_t endbin = x_finish - x_begin;
    const float *y_start  = y_begin + startbin;
    const float *y_finish = y_begin + endbin;
    
    m_group_counts.push_back( std::make_pair(y_start, y_finish) );
    m_group_energies.push_back( std::make_pair(x_start, x_finish) );
    
    m_nbins_used += (x_finish - x_start);
#if( WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL > 2 )
    cout << "Peak group " << i << " of " << m_grouped_candidates.size()
    << " has peaks with means { ";
    for( size_t j = 0; j < m_grouped_candidates[i].size(); ++j )
    cout << (j==0?"":", ") << m_grouped_candidates[i][j].mean();
    cout << " } and extends from " << lowere << " to " << uppere  << " keV"
    << endl;
#endif
  }//for( size_t i = 0; i < m_grouped_candidates.size(); ++i )
  
  
  //    cout << "Will fit for " << m_grouped_isfixed.size() << " groups of peaks."
  //         << endl;
  
  m_inited = true;
  return true;
}//bool init()


double AutoPeakSearchChi2Fcn::Up() const
{
  return 1.0;
}

double AutoPeakSearchChi2Fcn::operator()( const std::vector<double> &params ) const
{
  return eval_chi2( params );
}


double AutoPeakSearchChi2Fcn::peak_sigma( const double energy, const std::vector<double> &pars ) const
{
  if( pars.empty() )
    throw runtime_error( "peak_sigma: invalid intput" );
  
  double sigma = pars[0];
  
  if( m_resolution_type == SqrtEnergy )
  {
    if( pars.size() > 1 )
      sigma += (pars[1] * std::pow(energy, 0.5));
    
    if( pars.size() > 2 )
      sigma += pars[2] / (energy + 10.0);
  }else
  {
    for( size_t i = 1; i < pars.size(); ++i )
    sigma += (pars[i] * std::pow( energy, static_cast<double>(i) ));
  }
  
  const size_t channel = m_meas->find_gamma_channel(energy);
  sigma = max( sigma, (double)m_meas->gamma_channel_width(channel) );
  
  return sigma;
}//double peak_sigma(...)


ROOT::Minuit2::MnUserParameters AutoPeakSearchChi2Fcn::initial_parameters() const
{
  //Parameter order:
  //  Group 1, peak 1 mean
  //  Group 1, peak 2 mean
  //  Group 1, ...
  //  Group 1, peak N mean
  //  Group 2, peak 1 mean
  //  Group 2, peak 2 mean
  //  Group 2, ...
  //  Group 2, peak M mean
  //  ...
  //  Group K, peak J mean
  //  Peak Resolution parameter 1
  //  Peak Resolution parameter 2
  //  Peak Resolution parameter ...
  //  Peak Resolution parameter L
  
  if( !m_inited )
    throw runtime_error( "AutoPeakSearchChi2Fcn: you must call init before initial_parameters()" );
  
  ROOT::Minuit2::MnUserParameters pars;
  
  assert( m_grouped_candidates.size() == m_grouped_isfixed.size() );
  
  double maxsigma = 0.0;
  PeakDef leftpeak = m_grouped_candidates[0][0];
  PeakDef rightpeak = m_grouped_candidates[0][0];
  
  for( size_t group = 0; group < m_grouped_candidates.size(); ++group )
  {
    const vector<PeakDef> &peaks = m_grouped_candidates[group];
    const vector<bool> &isfixed = m_grouped_isfixed[group];
    
    assert( peaks.size() == isfixed.size() );
    
    for( size_t i = 0; i < peaks.size(); ++i )
    {
      char buffer[256];
      snprintf( buffer, sizeof(buffer), "mean_%i_%i", int(group), int(i) );
      
      const PeakDef &p = peaks[i];
      
      if( p.mean() < leftpeak.mean() )
        leftpeak = p;
      if( p.mean() > rightpeak.mean() )
        rightpeak = p;
      
      const double mean = p.mean();
      const double sigma = p.sigma();
      maxsigma = std::max( maxsigma, sigma );
      
      if( isfixed[i] )
        pars.Add( buffer, mean, 0.05*sigma, mean-0.2*sigma, mean+0.2*sigma );
      else
        pars.Add( buffer, mean, 0.2*sigma, mean-sigma, mean+sigma );
    }//for( size_t i = 0; i < peaks.size(); ++i )
  }//for( size_t group = 0; group < m_grouped_candidates.size(); ++group )
  
  size_t npeaks = 0;
  double meanwidth = 0.0;
  for( size_t group = 0; group < m_grouped_candidates.size(); ++group )
  {
    const vector<PeakDef> &peaks = m_grouped_candidates[group];
    for( size_t i = 0; i < peaks.size(); ++i )
    {
      ++npeaks;
      meanwidth += peaks[i].sigma();
    }//for( size_t i = 0; i < peaks.size(); ++i )
  }//for( size_t group = 0; group < m_grouped_candidates.size(); ++group )
  
  meanwidth /= npeaks;
  const double span = m_x->back() - m_x->front();
  const double meanbinwidth = span / m_x->size();
  
  switch( m_resolution_type )
  {
    case Polynomial0thOrder:
    {
      const double minwidth = m_isHPGe ? meanbinwidth : 4*meanbinwidth;
      pars.Add( "ResolutionZeroth", meanwidth, meanbinwidth,  minwidth, 2.0*maxsigma );
      
      break;
    }//case Polynomial0thOrder:
      
    case Polynomial1stOrder:
    case Polynomial2ndOrder:
    case Polynomial3rdOrder:
    {
      if( npeaks < 2 )
        throw std::logic_error( "Need more peaks for Polynomial1stOrder" );
      
      double de = rightpeak.mean() - leftpeak.mean();
      double slope = (rightpeak.sigma() - leftpeak.sigma()) / de;
      slope = std::max( slope, 0.0 );
      const double r0 = leftpeak.sigma() - slope*leftpeak.mean();
      pars.Add( "ResolutionZeroth", r0, meanbinwidth, 0.0, span/10.0 );
      pars.Add( "ResolutionFirst",  slope, meanbinwidth/span, 0.0, span/10.0 );
      
      if( m_resolution_type >= Polynomial2ndOrder )
        pars.Add( "ResolutionSecond", 0.0, meanbinwidth, -span/100.0, span/100.0 );
      
      if( m_resolution_type >= Polynomial2ndOrder )
        pars.Add( "ResolutionSecond", 0.0, meanbinwidth/span/span, -0.05/span/span, 0.05/span/span );
      
      if( m_resolution_type == Polynomial3rdOrder )
        pars.Add( "ResolutionThird", 0.0, meanbinwidth/span/span/span, -0.05/span/span/span, 0.05/span/span/span );
      
      break;
    }
      
    case SqrtEnergy:
    {
      if( npeaks < 2 )
        throw std::logic_error( "Need more peaks for SqrtEnergy" );
      
      //fwhm = a0 + a1*sqrt(e) + a3/(e+10)
      //fwhm(e_left) - a1*sqrt(e_left) = fwhm(e_right) - a1*sqrt(e_right)
      //(fwhm(e_left) - fwhm(e_right))/(sqrt(e_left) - sqrt(e_right)) = a2
      const double a1 = (leftpeak.sigma() - rightpeak.sigma()) / (sqrt(leftpeak.mean()) - sqrt(rightpeak.mean()));
      const double a0 = leftpeak.sigma() - a1*sqrt(leftpeak.sigma());
      
      pars.Add( "ResolutionZeroth", a1, meanbinwidth, 0.0, span/10.0 );
      pars.Add( "ResolutionFirst",  a0, meanbinwidth/span, -span/10.0, span/10.0 );
      
      break;
    }
      
  }//switch( m_resolution_type )
  
  
  return pars;
}


double AutoPeakSearchChi2Fcn::eval_chi2( const std::vector<double> &params ) const
{
  if( !m_inited )
    throw runtime_error( "AutoPeakSearchChi2Fcn: you must call init before eval_chi2()" );
  
  std::vector<PeakDef > peaks;
  
  return pars_to_peaks( peaks, params );
}//double eval_chi2( const std::vector<double> &params ) const


void AutoPeakSearchChi2Fcn::fit_peak_group( const vector<PeakDef> &peaks,
                                           const vector<double> &resolution_coefs,
                                           const size_t group,
                                           const double *pars,
                                           double &chi2,
                                           vector<PeakDef> &results ) const
{
  assert( peaks.size() );
  
  std::shared_ptr<const PeakContinuum> cont = peaks[0].continuum();
  
  const int num_polynomial_terms = ([&cont]() -> int {
    switch( cont->type() )
    {
      case PeakContinuum::NoOffset: case PeakContinuum::External:
        return 0;
        
      case PeakContinuum::Constant: case PeakContinuum::Linear:
      case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
        return cont->type() - PeakContinuum::NoOffset;
        
      case PeakContinuum::FlatStep:
      case PeakContinuum::LinearStep:
      case PeakContinuum::BiLinearStep:
        return 2 + (cont->type() - PeakContinuum::FlatStep);
    }//switch( cont->type() )

    assert(0);
    throw std::runtime_error( "Somehow invalid continuum polynomial type." );
    return 0;
  })();
  
  const bool isStepContinuum = ([&cont]() -> bool {
    switch( cont->type() )
    {
      case PeakContinuum::NoOffset: case PeakContinuum::External:
      case PeakContinuum::Constant: case PeakContinuum::Linear:
      case PeakContinuum::Quadratic: case PeakContinuum::Cubic:
        return false;
        
      case PeakContinuum::FlatStep:
      case PeakContinuum::LinearStep:
      case PeakContinuum::BiLinearStep:
        return true;
    }//switch( cont->type() )

    assert( 0 );
    throw std::runtime_error( "Somehow invalid continuum polynomial type." );
    return 0;
  })();
  
  
  std::vector<double> means, sigmas;
  
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    means.push_back( pars[i] );
    sigmas.push_back( peak_sigma( pars[i], resolution_coefs) );
  }
  
  const float *x_start = m_group_energies[group].first;
  const float *x_end = m_group_energies[group].second;
  const float *y_start = m_group_counts[group].first;
  
  const size_t nregionbin = x_end - x_start;
  
  const std::vector<PeakDef> fixedAmpPeaks;
  std::vector<double> amplitudes, amplitudes_uncerts;
  std::vector<double> continuum_coeffs, continuum_coeffs_uncerts;
  
#define USE_LLS_TO_FIT_MULTIPLE_PEAKS 1
  
#if( !USE_LLS_TO_FIT_MULTIPLE_PEAKS )
  if( means.size() > 1 )
  {
    double datasum = 0.0;
    for( const float *x = x_start; x != x_end; ++x )
    datasum += *x;
    
    const size_t lower_channel = m_meas->find_gamma_channel( *x_start );
    const size_t upper_channel = lowerbin + nregionbin;
    const int npeak = static_cast<int>( means.size() );
    const PeakContinuum::OffsetType type = cont->type();
    MultiPeakFitChi2Fcn chifcn( npeak, m_meas, type, lower_channel, upper_channel );
    
    
    ROOT::Minuit2::MnUserParameters inputPrams;
    if( type >= PeakContinuum::Constant )
      inputPrams.Add( "P0",  0.5*( (*x_start) + (*(x_end-1)) ) );
    if( type >= PeakContinuum::Linear )
      inputPrams.Add( "P1",  0.0 );
    if( type >= PeakContinuum::Quadratic )
      inputPrams.Add( "P2",  0.0 );
    if( type >= PeakContinuum::Cubic )
      inputPrams.Add( "P3",  0.0 );
    
    for( size_t i = 0; i < means.size(); ++i )
    {
      const string parnum = std::to_string(i);
      if( i == 0 )
        inputPrams.Add( "sigma" + parnum, sigmas[i] );
      else
        inputPrams.Add( "sigma" + parnum, -sigmas[i] );
      
      inputPrams.Add( "mean" + parnum, means[i] );
      inputPrams.Add( "amplitude" + parnum, 0.25*datasum, 0.1*datasum, 0, datasum );
    }//for( size_t i = 0; i < inpeaks.size(); ++i )
    
    ROOT::Minuit2::MnUserParameterState inputParamState( inputPrams );
    ROOT::Minuit2::MnStrategy strategy( 0 );
    ROOT::Minuit2::MnMinimize fitter( chifcn, inputParamState, strategy );
    ROOT::Minuit2::FunctionMinimum minimum = fitter( 2500, 1.0 );
    const vector<double> values = fitter.Parameters().Params();
    const vector<double> errors = fitter.Parameters().Errors();
    chi2 = chifcn.DoEval( &values[0], false );
    
    for( size_t i = 0; i < means.size(); ++i )
    {
      amplitudes.push_back( values[type+3*i+2] );
      amplitudes_uncerts.push_back( errors[type+3*i+2] );
    }
    
    continuum_coeffs.insert( continuum_coeffs.begin(), values.begin(), values.begin() + type );
    continuum_coeffs_uncerts.insert( continuum_coeffs_uncerts.begin(), errors.begin(), errors.begin() + type );
  }else
  {
#endif
    try
    {
      // TODO: currently this
      const PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;
      const double * const skew_parameters = nullptr;
      
      chi2 = fit_amp_and_offset( x_start, y_start, nregionbin,
                                num_polynomial_terms,
                                isStepContinuum,
                                cont->lowerEnergy(),
                                means, sigmas,
                                fixedAmpPeaks,
                                skew_type,
                                skew_parameters,
                                amplitudes,
                                continuum_coeffs,
                                amplitudes_uncerts,
                                continuum_coeffs_uncerts );
    }catch(...)
    {
      for( size_t i = 0; i < nregionbin; ++i )
      chi2 += y_start[i];
      amplitudes.clear();
      amplitudes.resize( means.size(), 0.0 );
      amplitudes_uncerts = amplitudes;
      continuum_coeffs = vector<double>( cont->type(), 0.0 );
      continuum_coeffs_uncerts = continuum_coeffs;
    }
#if( !USE_LLS_TO_FIT_MULTIPLE_PEAKS )
  }//if( means.size() > 1 ) / else
#endif
  
  std::shared_ptr<PeakContinuum> fitcont( new PeakContinuum() );
  fitcont->setRange( cont->lowerEnergy(), cont->upperEnergy() );
  fitcont->setType( cont->type() );
  fitcont->setParameters( cont->lowerEnergy(), continuum_coeffs,
                         continuum_coeffs_uncerts );
  for( size_t i = 0; i < amplitudes.size(); ++i )
  {
    PeakDef peak;
    peak.setContinuum( fitcont );
    peak.setMean( means[i] );
    peak.setAmplitude( std::max(amplitudes[i],0.0) );
    peak.setSigma( sigmas[i] );
    peak.setAmplitudeUncert( amplitudes_uncerts[i] );
    
    results.push_back( peak );
  }//for( size_t i = 0; i < amplitudes.size(); ++i )
}//fit_peak_group()


//pars_to_peaks(): returns chi2
double AutoPeakSearchChi2Fcn::pars_to_peaks( std::vector<PeakDef> &resultpeaks,
                                            const std::vector<double> &pars,
                                            const std::vector<double> &errors ) const
{
  resultpeaks.clear();
  
  if( !m_inited )
    throw runtime_error( "AutoPeakSearchChi2Fcn: you must call init before pars_to_peaks()" );
  
  if( m_resolution_type != SqrtEnergy )
    assert( pars.size() == (m_num_peakstotal+m_resolution_type+1) );
  else
    assert( pars.size() == (m_num_peakstotal+2) );
  
  vector<double> rescoef( pars.begin()+m_num_peakstotal, pars.end() );
  
  const size_t ngroups = m_grouped_candidates.size();
  vector<double> chi2s( ngroups, 0.0 );
  vector< vector<PeakDef> > peakdefs( ngroups );
  
  size_t peakn = 0;
  
  SpecUtilsAsync::ThreadPool pool;
  
  for( size_t group = 0; group < m_grouped_candidates.size(); ++group )
  {
    const vector<PeakDef> &peaks = m_grouped_candidates[group];
    
    const double *startpars = &(pars[0]) + peakn;
    
    pool.post( boost::bind( &AutoPeakSearchChi2Fcn::fit_peak_group, this,
                           boost::cref(peaks), boost::cref(rescoef), group, startpars,
                           boost::ref(chi2s[group]), boost::ref(peakdefs[group]) ) );
    
    peakn += peaks.size();
  }//for( size_t group = 0; group < m_grouped_candidates.size(); ++group )
  
  pool.join();
  
  double chi2 = 0.0;
  for( size_t i = 0; i < chi2s.size(); ++i )
  {
    chi2 += chi2s[i];
    for( size_t j = 0; j < peakdefs[i].size(); ++j )
    resultpeaks.push_back( peakdefs[i][j] );
  }
  
  std::sort( resultpeaks.begin(), resultpeaks.end(), &PeakDef::lessThanByMean );
  chi2 += punishment( resultpeaks );
  
  return chi2;
}//double pars_to_peaks(...)


//punishment(..): the punishment to the chi2 for peaks being to close together
//  or not statistically significant
double AutoPeakSearchChi2Fcn::punishment( const std::vector<PeakDef> &peaks ) const
{
  vector<double> means( peaks.size() ), sigmas( peaks.size() ), amps( peaks.size() );
  for( size_t i = 0; i < peaks.size(); ++i )
  {
    means[i] = peaks[i].mean();
    sigmas[i] = peaks[i].sigma();
    amps[i] = peaks[i].amplitude();
  }//for( size_t i = 0; i < m_npeak; ++i )
  
  const double punishment_chi2 = 2.0*m_nbins_used;
  
  double chi2 = closeness_punishment( means, sigmas );
  
  //If the peak area is statistically insignificant on the interval
  //  -1.75sigma to 1.75, then punish!
  for( size_t i = 0; i < means.size(); ++i )
  {
    const double lower_energy = means[i] - 1.75*sigmas[i];
    const double upper_energy = means[i] + 1.75*sigmas[i];
    
    const size_t binstart = lower_bound( m_x->begin(), m_x->end(), lower_energy )
    - m_x->begin();
    size_t binend = lower_bound( m_x->begin(), m_x->end(), upper_energy )
    - m_x->begin();
    binend = std::min( binend, m_y->size() );
    
    double dataarea = 0.0;
    for( size_t bin = binstart; bin < binend; ++bin )
    dataarea += (*m_y)[bin];
    
    if( peaks[i].amplitude() < 2.0*sqrt(dataarea) )
    {
      if( amps[i] <= 1 )
        chi2 += 100.0 * punishment_chi2;
      else
        chi2 += 0.5*(sqrt(dataarea)/amps[i]) * punishment_chi2;
    }//if( peaks[i].amplitude() < 2.0*sqrt(dataarea) )
  }//for( size_t i = 0; i < peaks.size(); ++i )
  
  return chi2;
}//punishment(...)


double AutoPeakSearchChi2Fcn::closeness_punishment( const vector<double> &means,
                                                   const vector<double> &sigmas ) const
{
  double chi2 = 0.0;
  const double punishment_chi2 = 2.0*m_nbins_used;
  
  //punish if the peaks are too close.
  for( size_t i = 1; i < means.size(); ++i )
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
}//closeness_punishment(...)

//other_peaks may contain peak
bool AutoPeakSearchChi2Fcn::significance_test( const PeakDef &peak,
                                              std::vector<PeakDef> other_peaks,
                                              double *chi2DOF,
                                              double *grosCountsSigma
                                              ) const
{
  if( !peak.gausPeak() )
    throw runtime_error( "significance_test can only evaluate gaussian defined peaks" );
  
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  for( size_t i = 0; i < other_peaks.size(); ++i )
  if( continuum != other_peaks[i].continuum() )
    throw runtime_error( "significance_test all input peaks must share a continuum" );
  
  std::vector<PeakDef>::iterator pos = other_peaks.end();
  for( size_t i = 0; i < other_peaks.size(); ++i )
  {
    const PeakDef &other = other_peaks[i];
    
    if( fabs(other.mean() - peak.mean()) < 0.0001
       && fabs(other.sigma() - peak.sigma()) < 0.0001 )
    {
      pos = other_peaks.begin() + i;
      break;
    }
  }//for( size_t i = 0; i < other_peaks.size(); ++i )
  
  //  pos = std::find( other_peaks.begin(), other_peaks.end(), peak );
  if( pos != other_peaks.end() )
    other_peaks.erase( pos );
  
  //Some basic quality checks
  if( IsNan(peak.mean()) || IsInf(peak.mean()) )
    return false;
  if( IsNan(peak.sigma()) || IsInf(peak.sigma()) )
    return false;
  if( IsNan(peak.amplitude()) || IsInf(peak.amplitude()) )
    return false;
  if( peak.amplitude() <= 0.0 || peak.sigma() <= 0.0 )
    return false;
  
  double withChi2 = 0.0, withoutChi2 = 0.0;
  
  const float *xbegin = &((*m_x)[0]);
  const float *xend = xbegin + m_x->size();
  
  const double lowerEnergy = peak.continuum()->lowerEnergy();
  const double upperEnergy = peak.continuum()->upperEnergy();
  
  const float *begin = std::upper_bound( xbegin, xend, lowerEnergy );
  if( begin > xbegin )
    --begin;
  const float *end = std::upper_bound( xbegin, xend, upperEnergy );
  if( end < xend )
    ++end;
  if( end > (xend-1) )
    end = (xend-1);
  
  if( end <= begin )
    throw std::logic_error( "significance_test: error finding range" );
  
  const size_t endbin = end - xbegin;
  const size_t startbin = begin - xbegin;
  for( size_t bin = startbin; bin < endbin; ++bin )
  {
    const float loweredge = (*m_x)[bin];
    const float upperedge = (*m_x)[bin+1];
    const double continuumarea = continuum->offset_integral( loweredge, upperedge, m_meas );
    const double peakarea = peak.gauss_integral( loweredge, upperedge );
    double otherPeakArea = 0.0;
    for( size_t i = 0; i < other_peaks.size(); ++i )
    otherPeakArea += other_peaks[i].gauss_integral( loweredge, upperedge );
    const double dataarea = std::max( (*m_y)[bin], 1.0f );
    
    withChi2 += std::pow(peakarea+otherPeakArea+continuumarea-dataarea, 2.0) / dataarea;
    withoutChi2 += std::pow(otherPeakArea+continuumarea-dataarea, 2.0) / dataarea;
  }//for( size_t bin = startbin; bin < endbin; ++bin )
  
  withChi2 /= (endbin-startbin);
  withoutChi2 /= (endbin-startbin);
  
  const double deltaChi2 = withoutChi2 - withChi2;
  
  
  
  begin = std::upper_bound( xbegin, xend, peak.mean()-1.18*peak.sigma() );
  if( begin > xbegin )
    --begin;
  end = std::upper_bound( xbegin, xend, peak.mean()+1.18*peak.sigma() );
  if( end < xend )
    ++end;
  if( end > (xend-1) )
    end = (xend-1);
  
  const size_t beingint = begin-xbegin;
  const size_t endint = end-xbegin;
  
  double datamidarea = 0.0;
  for( size_t i = beingint; i < endint; ++i )
  datamidarea += (*m_y)[i];
  
  const double gausmidarea = peak.gauss_integral( *begin, *end );
  
  const double gross_counts_sig = gausmidarea/sqrt(datamidarea);
  
  if( chi2DOF )
    *chi2DOF = deltaChi2;
  if( grosCountsSigma )
    *grosCountsSigma = gross_counts_sig;
  
  //    double *chi2DOF = 0,
  //    double *grosCountsSigma = 0
  
  //    cerr << "Peak at mean=" << peak.mean() << ", sigma=" << peak.sigma() << " has deltaChi2=" << deltaChi2
  //         << ", gausmidarea=" << gausmidarea << ", datamidarea=" << datamidarea
  //         << ", (gausmidarea/sqrt(datamidarea))=" << (gausmidarea/sqrt(datamidarea))
  //         << endl;
  
  return ((deltaChi2 > m_min_chi2_dof_thresh)
          && gross_counts_sig>=m_min_gross_counts_sig_thresh);
}//bool significance_test( ... )



std::vector<PeakDef> AutoPeakSearchChi2Fcn::filter_peaks( const std::vector<PeakDef> &all_peaks ) const
{
  std::vector<PeakDef> keepers;
  
  typedef std::map< std::shared_ptr<const PeakContinuum>, vector<PeakDef> > peakmap_t;
  peakmap_t peakmap;
  
  for( size_t i = 0; i < all_peaks.size(); ++i )
  peakmap[all_peaks[i].continuum()].push_back( all_peaks[i] );
  
  for( const peakmap_t::value_type &t : peakmap )
  {
    const vector<PeakDef> &v = t.second;
    for( size_t i = 0; i < v.size(); ++i )
    if( significance_test(v[i], v) )
      keepers.push_back( v[i] );
  }
  
  std::sort( keepers.begin(), keepers.end(), &PeakDef::lessThanByMean );
  
  return keepers;
}//filter_peaks(...)


bool find_spectroscopic_extent( std::shared_ptr<const Measurement> meas,
                               size_t &lower_channel,
                               size_t &upper_channel )
{
  if( !meas )
    return false;
  
  const vector<float> &channel_counts = *meas->gamma_counts();
  const size_t nbin = channel_counts.size();
  
  
  //First detect where spectrum begins
  const int side_bins = 5;
  const int order = 2;
  const int derivative = 2;
  vector<float> smoothed_2nd;
  smoothSpectrum( meas, side_bins, order, derivative, smoothed_2nd );
  
  size_t channel = 0;
  while( channel < nbin && smoothed_2nd[channel]>-1.0 )
    ++channel;
  
  //    const size_t first_negative_2nd_deriv = channel;
  
  while( channel < nbin && smoothed_2nd[channel]<0.05 )
    ++channel;
  //  cerr << "Lower channel = " << meas->channel_energies()->at(channel) << endl;
  lower_channel = std::min(channel-0,smoothed_2nd.size()-1);
  
  if( lower_channel > (nbin/3) )
  {
    //cout << "The lower threshold bin (" << lower_channel << " of " << nbin
    //<< ") is to high, skipping file for further analysis." << endl;
    lower_channel = upper_channel = 0;
    return false;
  }//
  
  size_t upperlastchannel = nbin - 1;
  while( upperlastchannel > 0 && channel_counts[upperlastchannel] < 5.0f )
    --upperlastchannel;
  
  //Start at the turn on energy, but make sure were in a region of the spectra
  //  with decent statistics (e.g. at least 20 counts per bin)
  size_t lastchannel = lower_channel;
  while( lastchannel < nbin && channel_counts[lastchannel] < 20.0 )
    ++lastchannel;
  
  //Now find where the spectrum drops below X counts per channel for a number
  //  of channels in a row
  int numbelow = 0;
  const int nbelowlimit = std::max( int(std::ceil(0.0015*nbin)), 3 );
  while( lastchannel < (nbin-1) && numbelow <= nbelowlimit )
    numbelow = (channel_counts[++lastchannel] > 0.1) ? 0 : numbelow + 1;
  
  upper_channel = (lastchannel == nbin) ? size_t(nbin-1) : lastchannel;
  upper_channel = std::max( upper_channel, upperlastchannel+1 );
  upper_channel = std::min( upper_channel, nbin-1 );
  
  return true;
}//bool find_spectroscopic_extent(...)

std::vector<PeakDef> search_for_peaks( const std::shared_ptr<const Measurement> meas,
                                      const std::vector<PeakDef> &origpeaks,
                                      const bool isHPGe )
{
  if( !meas || !meas->gamma_counts() )
    return origpeaks;
  
  const double min_chi2_dof_thresh         = isHPGe ? 0.2   : 3.5;
  const double min_gross_counts_sig_thresh = isHPGe ? 3.0   : 3;
  const double above_line_chi2_thresh      = isHPGe ? 4.075 : 4.075;
  const int side_bins                      = isHPGe ? 7     : 10;
  const int smooth_order                   = isHPGe ? 3     : 2;
  const double second_deriv_thresh         = isHPGe ? -0.02 : 0.04;
  const double stat_thresh                 = isHPGe ? 1.0   : 2;
  const double width_thresh                = isHPGe ? 0.0   : 0.0;
  
  
  //initial_above_line_chi2_thresh: 3.5,
  //initial_second_deriv_thresh: 0.05,
  //initial_stat_thresh: 1.3,
  //final_min_chi2_dof_thresh: 1.1,
  //final_min_gross_counts_sig_thresh: 4.5
  
  return search_for_peaks( meas, min_chi2_dof_thresh,
                          min_gross_counts_sig_thresh,
                          above_line_chi2_thresh, side_bins, smooth_order,
                          second_deriv_thresh, stat_thresh, width_thresh,
                          origpeaks, isHPGe
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
                          , std::shared_ptr<const DetectorPeakResponse>()
#endif
                          );
  
}//search_for_peaks

std::vector<PeakDef> search_for_peaks( const std::shared_ptr<const Measurement> meas,
                                      const double min_chi2_dof_thresh,
                                      const double min_gross_counts_sig_thresh,
                                      const double above_line_chi2_thresh,
                                      const int side_bins,
                                      const int smooth_order,
                                      const double second_deriv_thresh,
                                      const double stat_thresh,
                                      const double width_thresh,
                                      const std::vector<PeakDef> &origpeaks, /*included in result, unmodified, wont have duplciate */
                                      const bool isHPGe
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
                                      , std::shared_ptr<const DetectorPeakResponse> detector
#endif
)
{
  
  vector<PeakDef> finalpeaks;
  
  if( !meas )
    throw runtime_error( "search_for_peaks: invalid input" );
  
  // TODO: currently doesn't account for/fit peak skew
  
  AutoPeakSearchChi2Fcn chi2fcn( meas, origpeaks, isHPGe );
  chi2fcn.m_min_chi2_dof_thresh = min_chi2_dof_thresh;
  chi2fcn.m_min_gross_counts_sig_thresh = min_gross_counts_sig_thresh;
  chi2fcn.m_nsigma_near_group = 3.0;
  
  
  chi2fcn.m_side_bins = side_bins;
  chi2fcn.m_smooth_order = smooth_order;
  chi2fcn.m_second_deriv_thresh = second_deriv_thresh;
  chi2fcn.m_stat_thresh = stat_thresh;
  chi2fcn.m_width_thresh = width_thresh;
  
  const bool success = chi2fcn.init();  //this is where the initial candidate peaks are actually produced
  
  if( !success )
    throw runtime_error( "Unable to init peak searching" );
  
  assert( success );
  
  
  ROOT::Minuit2::MnUserParameters params = chi2fcn.initial_parameters();
  ROOT::Minuit2::MnUserParameterState inputParamState( params );
  ROOT::Minuit2::MnStrategy strategy( 2 );
  
  ROOT::Minuit2::CombinedMinimizer fitter;
  unsigned int maxFcnCall = 0;
  const double tolerance = 0.001;
  
  
  ROOT::Minuit2::FunctionMinimum minimum
  = fitter.Minimize( chi2fcn, params, strategy, maxFcnCall, tolerance );
  
  params = minimum.UserState().Parameters();
  const vector<double> pars = params.Params();
  const vector<double> errors = params.Errors();
  
  std::vector<PeakDef> resultpeaks;
  
#if( WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL > 1 )
  const double fitchi2 = chi2fcn.pars_to_peaks( resultpeaks, pars, errors );
  cerr << "Fit " << resultpeaks.size() << " peaks w/ chi2=" << fitchi2 << endl;
  for( size_t i = 0; i < resultpeaks.size(); ++i )
  cerr << "\tunfiltered: mean=" << resultpeaks[i].mean() << ", area=" << resultpeaks[i].amplitude() << endl;
#endif
  
  std::vector<PeakDef> filtered = chi2fcn.filter_peaks( resultpeaks );
  
#if( WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL > 1 )
  cerr << "After Filtering there are still " << filtered.size() << " peaks" << endl;
  for( size_t i = 0; i < filtered.size(); ++i )
  {
    cerr << "\tFiltered: mean=" << filtered[i].mean() << ", area=" << filtered[i].amplitude() << endl;
  }
#endif
  
  //Here we should fit for the FWHM funcitonal form if we have a ny peaks
  finalpeaks = filtered;
  vector<double> rescoef( pars.begin()+chi2fcn.num_initial_candidate_peaks(), pars.end() );
  
  if( filtered.size() )
  {
    //get the FWHM as a function of energy
    //      chi2fcn.m_min_gross_counts_sig_thresh *= 2.0;
    //      chi2fcn.m_min_chi2_dof_thresh *= 2.0;
    
    //Scan the spectrum over the valid range to find candidate peaks, filtering as we go
    //  -could do something really simple like have a linear continuum that touches
    //   at +-2 sigma; If there is another peak there, than have it touch its continuum.
    //   If the
    vector<double> chi2aboveline, chi2abovelineEnergies, chi2abovelineareas;
    vector< vector<double> > chi2abovelinecontcoefs;
    vector< double > chi2abovelinecontstart;
    
    const vector<float> &energies = *meas->channel_energies();
    const vector<float> &spectrum = *meas->gamma_counts();
    
    const float lowerenergy = energies.at( chi2fcn.lower_spectrum_channel() );
    const float upperenergy = energies.at( chi2fcn.upper_spectrum_channel() );
    
#if( WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL > 1 )
    cerr << "Will scan from " << lowerenergy << " to " << upperenergy << " keV" << endl;
    cerr << "Inital sigma=" << chi2fcn.peak_sigma(lowerenergy,rescoef) << endl;
    cerr << "Upper sigma=" << chi2fcn.peak_sigma(upperenergy,rescoef) << endl;
#endif
    
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
    ofstream spectrumfile("spectrum.csv");
    ofstream chi2peakfile("chi2withpeak.csv");
    ofstream chi2nopeakfile("chi2nopeak.csv");
    ofstream chi2ratiofile("chi2ratio.csv");
    ofstream chi2difffile("chi2diff.csv");
    ofstream linefitfile("linefit.csv");
    ofstream peakslinefitfile("peakslinefit.csv");
    ofstream areaabovelinefile("areaaboveline.csv");
    ofstream issigfile("issig.csv");
    
    for( size_t i = 0; i < spectrum.size(); ++i )
    spectrumfile << energies[i] << "," << spectrum[i] << endl;
    
    if( !!detector && detector->hasResolutionInfo() )
    {
      ofstream foundfwhm("foundFWHM.csv");
      ofstream gadrasfwhm("gadrasFWHM.csv");
      for( size_t i = 0; i < spectrum.size(); ++i )
      {
        foundfwhm << energies[i] << "," << chi2fcn.peak_sigma(energies[i],rescoef) << endl;
        gadrasfwhm << energies[i] << "," << detector->peakResolutionSigma(energies[i]) << endl;
      }
    }
#endif
    
    
    float deltae = chi2fcn.peak_sigma(lowerenergy,rescoef);
    for( float energy = lowerenergy + deltae;
        energy < upperenergy - deltae;
        energy += deltae )
    {
      
      deltae = 0.25*chi2fcn.peak_sigma(energy,rescoef);
      const size_t channel = meas->find_gamma_channel(energy);
      const float binwidth = meas->gamma_channel_width(channel);
      if( deltae < binwidth )
        deltae = 1.5*binwidth;
      
      const double sigma = chi2fcn.peak_sigma(energy,rescoef);
      
      //Fit for a peak
      vector<PeakDef> existingpeaks = peaksTouchingRange( energy-1.5*sigma, energy+1.5*sigma, finalpeaks );
      
      double lower_range = energy - 2.0*sigma;
      double upper_range = energy + 2.0*sigma;
      
      bool isInOtherPeak = false;
      
      for( const PeakDef &p : existingpeaks )
      {
        lower_range = std::min( lower_range, p.lowerX() );
        upper_range = std::max( upper_range, p.upperX() );
        
        isInOtherPeak = (isInOtherPeak || (fabs(p.mean()-energy) < 2.0*p.sigma()));
      }
      
      if( isInOtherPeak )
      {
        //          cout << energy << " is in other peak, but not skipping" << endl;
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
        chi2peakfile << energy << "," << 0.2 << endl;
        chi2nopeakfile << energy << "," << 0.2 << endl;
        
        issigfile        << energy << "," << 0 << endl;
        chi2peakfile     << energy << "," << 0 << endl;
        chi2nopeakfile   << energy << "," << 0 << endl;
        peakslinefitfile << energy << "," << 0 << endl;
        linefitfile      << energy << "," << 0 << endl;
#endif
        continue;
      }
      
      try
      {
        const std::vector<PeakDef> fixedAmpPeaks;
        vector<double> amplitudes, continuum_coeffs, amplitudes_uncerts, continuum_coeffs_uncerts;
        
        const size_t startbin = lower_bound( energies.begin(), energies.end(), lower_range ) - energies.begin();
        const size_t endbin = lower_bound( energies.begin(), energies.end(), upper_range ) - energies.begin();
        
        assert( startbin != energies.size() );
        assert( endbin >= startbin );
        
        const float *x_start = &(energies[0]) + startbin;
        const float *data = &(spectrum[0]) + startbin;
        const size_t nbin = endbin - startbin;
        
        vector<double> means, sigmas;
        for( const PeakDef &p : existingpeaks )
        {
          //            cout << "Not adding to existing peaks" << endl;
          means.push_back( p.mean() );
          sigmas.push_back( p.sigma() );
        }
        
        means.push_back( energy );
        sigmas.push_back( sigma );
        
        int num_polynomial_terms = 2;  //linear
        const bool step_continuum = false;
        /*
         if( means.size() > 2 )
         num_polynomial_terms = 3;
         if( nbin > 15 )  //15 is arbitrarily chosen right now
         num_polynomial_terms = 3;
         */
        
        // TODO: currently doesn't account for/fit peak skew
        const PeakDef::SkewType skew_type = PeakDef::SkewType::NoSkew;
        const double * const skew_pars = nullptr;
        
        const double chi2 =  fit_amp_and_offset( x_start, data, nbin,
                                                num_polynomial_terms,
                                                step_continuum,
                                                energy, means, sigmas,
                                                fixedAmpPeaks,
                                                skew_type,
                                                skew_pars,
                                                amplitudes, continuum_coeffs,
                                                amplitudes_uncerts, continuum_coeffs_uncerts );
        
        //we have amplitudes, uncertainties
        const double fitamp = amplitudes.back();
        const double fitampuncert = amplitudes_uncerts.back();
        
        double withpeakchi2 = 0.0, withoutpeakchi2 = 0.0;
        for( size_t bin = 0; bin < nbin; ++bin )
        {
          const double x0 = x_start[bin];
          const double x1 = x_start[bin+1];
          double y_pred = 0.0;
          for( size_t col = 0; col < continuum_coeffs.size(); ++col )
          {
            const double exp = col + 1.0;
            const double x0cont = x0 - energy;
            const double x1cont = x1 - energy;
            y_pred += continuum_coeffs[col] * (1.0/exp) * (pow(x1cont,exp) - pow(x0cont,exp));
          }//for( int order = 0; order < maxorder; ++order )
          
          const double uncert = (data[bin] > MIN_CHANNEL_UNCERT ? sqrt( data[bin] ) : 1.0);
          
          if( y_pred < 0.0 )
            y_pred = 0.0;
          
          for( size_t i = 0; i < (amplitudes.size()-1); ++i )
          y_pred += amplitudes[i]*PeakDists::gaussian_integral( means[i], sigmas[i], x0, x1 );
          
          for( size_t i = 0; i < fixedAmpPeaks.size(); ++i )
          y_pred += fixedAmpPeaks[i].gauss_integral( x0, x1 );
          
          withoutpeakchi2 += std::pow( (y_pred - data[bin]) / uncert, 2.0 );
          
          y_pred += amplitudes.back()*PeakDists::gaussian_integral( means.back(), sigmas.back(), x0, x1 );
          withpeakchi2 += std::pow( (y_pred - data[bin]) / uncert, 2.0 );
        }
        
        double withpeakchi2dof = withpeakchi2 / (nbin+3);
        double withoutpeakchi2dof = withoutpeakchi2 / (nbin+2);
        
        if( fabs(withpeakchi2-chi2) < (0.05*chi2) || fabs(withpeakchi2-chi2)<0.01 )
        {
          //            cerr << "\n\nFailed chi2 check: withpeakchi2=" << withpeakchi2
          //                 << ", chi2=" << chi2 << endl;
          continue;
        }
        
        //          const double dataarea = gamma_integral( meas, energy-2.0*sigma, energy+2.0*sigma );
        
        double linechi2dof = 0.0;
        double peakscontheight = 0.0, lineonlyheight = 0.0;
        
        {
          //Here we could try fitting to a straight line, and try to use this to
          //  threshold off of as well....
          vector<float> lineonlyx, fakedata;
          for( size_t i = 0; i < nbin; ++i )
          lineonlyx.push_back( x_start[i] - energy );
          
          std::vector<double> lineonlycoeffs, lineonlyuncerts;
          const double linechi2 = fit_to_polynomial( &lineonlyx[0], data,
                                                    lineonlyx.size(), 1,
                                                    lineonlycoeffs,
                                                    lineonlyuncerts );
          linechi2dof = linechi2/(lineonlyx.size()-2);
          lineonlyheight = lineonlycoeffs[0];
          
          size_t midenergybin = lower_bound( lineonlyx.begin(), lineonlyx.end(), 0.0 ) - lineonlyx.begin();
          for( size_t col = 0; col < continuum_coeffs.size(); ++col )
          {
            const double exp = col + 1.0;
            const double x0cont = lineonlyx[midenergybin];
            const double x1cont = lineonlyx[midenergybin+1];
            peakscontheight += continuum_coeffs[col] * (1.0/exp) * (pow(x1cont,exp) - pow(x0cont,exp));
          }//for( int order = 0; order < maxorder; ++order )
        }
        
        
        PeakDef peak;
        peak.setMean( energy );
        peak.setSigma( sigma );
        peak.setAmplitude( fitamp );
        peak.continuum()->setType( PeakContinuum::Linear );
        peak.continuum()->setParameters( energy, &continuum_coeffs[0], &continuum_coeffs_uncerts[0] );
        double testChi2DOF, testGrosCountsSigma;
        vector<PeakDef> neighbors;
        /*const bool sig = */chi2fcn.significance_test( peak, neighbors,
                                                       &testChi2DOF, &testGrosCountsSigma );
        
        
        if( ((fitamp/fitampuncert) < 2.5)
           || (fitamp < 0.0)
           || (withoutpeakchi2dof-withpeakchi2dof) < 1.0
           || withpeakchi2dof > 2.0
           || withoutpeakchi2dof < 1.25
           /*|| (fitamp/sqrt(dataarea)<4)*/ )
        {
          withpeakchi2dof = 0.0;
          linechi2dof = 0.0;
          peakscontheight = 0.0;
          lineonlyheight = 0.0;
        }else
        {
          //fuck it, lets add this peak to
          finalpeaks.push_back( peak );
          
          /*
           const double sig = fitamp / fitampuncert;
           chi2aboveline.push_back( sig );
           chi2abovelineEnergies.push_back( energy );
           chi2abovelineareas.push_back( fitamp );
           chi2abovelinecontcoefs.push_back( vector<double>() );
           */
          //Then filter them, and refit peaks.  Should we
          
          //refitPeaksThatShareROI(...);
          
          //Write a new chi2 class that takes all the candidate peaks,
          //  and then decides which peaks to group together, and what
          //  detector response function to to use, then it fits for all
          //  the peaks using groupings.  Minuit fit parameters will be
          //  detector resolution function and peak energies.
          //  Can even waite to decide wich peaks are grouped together until
          //  parameters are specified.
          //  Once the new peaks are fuly fit for, then can filter them an
          //  additional time.
          //  Once thats all implemented, then can optimize all the
          //  parameters.
        }
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
        issigfile        << energy << "," << fitamp << endl;
        chi2peakfile     << energy << "," << withpeakchi2dof << endl;
        chi2nopeakfile   << energy << "," << linechi2dof     << endl;
        peakslinefitfile << energy << "," << peakscontheight << endl;
        linefitfile      << energy << "," << lineonlyheight  << endl;
#endif
        
        if( false )
        {//begin computation of straight line, skipping the peak area
          //next thing to try, fit a straight line, using about 1 sigma of data,
          //  then ~3.5 sigma not in the fit, then another 1 sigma of data, and
          //  look for the excess in the middle part, and use this to see if a
          //  peak is there (for HPGe, might need something else for NaI).
          vector<float> fakex, fakedata;
          const double fake_lower_range = energy - 2.5*sigma;
          const double fake_upper_range = energy + 2.5*sigma;
          
          const size_t fakestartbin = lower_bound( energies.begin(), energies.end(), fake_lower_range ) - energies.begin();
          const size_t fakeendbin = lower_bound( energies.begin(), energies.end(), fake_upper_range-sigma ) - energies.begin();
          
          size_t lastlowerbin = 0;
          for( size_t i = fakestartbin; fakex.empty() || energies[i] <= (fake_lower_range+sigma); ++i )
          {
            fakex.push_back( energies[i] );
            fakedata.push_back( spectrum[i] );
            lastlowerbin = i;
          }
          
          for( size_t i = fakeendbin; i < (fakeendbin+1) || energies[i] <= fake_upper_range; ++i )
          {
            fakex.push_back( energies[i] );
            fakedata.push_back( spectrum[i] );
          }
          
          const double fakeenergyoffset = fakex[0];
          chi2abovelinecontstart.push_back( fakeenergyoffset );
          for( size_t i = 0; i < fakex.size(); ++i )
          fakex[i] -= fakeenergyoffset;
          
          std::vector<double> fakepoly_coeffs, fakecoeff_uncerts;
          /*const double nomiddlefit = */ fit_to_polynomial( &fakex[0], &fakedata[0], fakex.size(), 1,
                                                            fakepoly_coeffs,
                                                            fakecoeff_uncerts );
          
          //now integrate the data for energy+-1.5 sigma
          //            const double middataarea = gamma_integral( meas, energy-1.5*sigma, energy+1.5*sigma );
          
          double midcontinuumarea = 0.0, middataarea = 0.0, midcontinuumareaplus = 0.0;
          const size_t midnbins = fakeendbin - lastlowerbin;
          for( size_t bin = lastlowerbin + 0.2*midnbins;
              bin < (fakeendbin-0.2*midnbins); ++bin )
          {
            double y_pred = 0.0, y_pred_plus = 0.0;
            for( size_t i = 0; i < fakepoly_coeffs.size(); ++i )
            {
              y_pred += fakepoly_coeffs[i] * std::pow( double(energies[bin]-fakeenergyoffset), double(i) );
              //                y_pred_plus += (fakepoly_coeffs[i]+fakecoeff_uncerts[i]) * std::pow( double(energies[bin]-fakeenergyoffset), double(i) );
            }
            
            y_pred = std::max( y_pred, 0.0 );
            y_pred_plus = y_pred + fakepoly_coeffs[0];
            //              y_pred_plus = std::max( y_pred_plus, 0.0 );
            
            midcontinuumarea += y_pred;
            midcontinuumareaplus += y_pred_plus;
            middataarea += spectrum[bin];
          }
          
          //            const double midcontinuumarea = (3.0*fakepoly_coeffs[0] + 0.5*(fakepoly_coeffs[1]*( pow(energy+1.5*sigma,2.0))-pow(energy-1.5*sigma,2.0))) / (3.0*sigma);
          
          //            cout << "nomiddlefit=" << nomiddlefit
          //                 << ", for " << fakex.size() << " bins" << endl;
          //            cout << "energy=" << energy
          //                 << ", midcontinuumarea=" << midcontinuumarea
          //                 << ", middataarea=" << middataarea << endl;
          //            const double continuum_uncert = fabs(midcontinuumareaplus - midcontinuumarea);
          
          
          const double midpeakarea = std::max(middataarea-midcontinuumarea,0.0);
          const double uncert = sqrt(middataarea); //sqrt(continuum_uncert*continuum_uncert + middataarea);
          
          const double sig = midpeakarea / uncert;
          
          //double sigmamidarea = std::max(middataarea-midcontinuumarea,0.0) / std::max( sqrt(middataarea), 1.0);
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
          areaabovelinefile << energy << "," << sig << endl;
#endif
          
          //            if( energy > 2102.5 && energy < 2103. )
          //            {
          //              cout << "found it sig=" << sig << endl;
          //            }
          
          chi2aboveline.push_back( sig /*sigmamidarea*/ );
          chi2abovelineEnergies.push_back( energy );
          chi2abovelineareas.push_back( middataarea-midcontinuumarea );
          chi2abovelinecontcoefs.push_back( fakepoly_coeffs );
        }//end computation of straight line, skipping the peak area
        
        
      }catch( std::exception &e )
      {
        cout << "Caught exception at " << energy << ": " << e.what() << endl;
        
      }
    }//for( loop over enrgies )
    
#if( WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL > 0 )
    cout << "chi2aboveline.size=" << chi2aboveline.size() << endl;
#endif
    
    while( 1 )
    {
      vector<double>::iterator maxeleiter = std::max_element( chi2aboveline.begin(), chi2aboveline.end() );
      
      if( maxeleiter==chi2aboveline.end() || (*maxeleiter) < above_line_chi2_thresh )
        break;
      
      const size_t index = maxeleiter - chi2aboveline.begin();
      //const double maxval = *maxeleiter;
      
      const double energy = chi2abovelineEnergies[index];
      const double area = chi2aboveline[index];
      
#if( WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL > 1 )
      cerr << "Adding peak to finalpeaks at energy=" << energy
      << ", with area " << area << " and sigma value " << maxval << endl;
#endif
      
      const double sigma = chi2fcn.peak_sigma(energy,rescoef);
      
      double wsum = 0.0;
      //zero out +-2 sigma of this max value
      for( size_t i = index; chi2abovelineEnergies[i] > (energy-2.0*sigma) && i > 1; --i )
      {
        wsum += chi2aboveline[i];
        chi2aboveline[i] = 0.0;
      }
      
      for( size_t i = index; chi2abovelineEnergies[i] < (energy+2.0*sigma) && i < chi2abovelineEnergies.size(); ++i )
      {
        wsum += chi2aboveline[i];
        chi2aboveline[i] = 0.0;
      }
      
#if( WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL > 1 )
      cout << "\twsum=" << wsum << ", area=" << area << endl;
#endif
      
      if( wsum < 4.0*area )
        continue;
      
      
      //Now go through and fit this peak, and add it to finalpeaks...
      PeakDef newpeak;
      newpeak.setMean( energy );
      newpeak.setSigma( sigma );
      newpeak.setAmplitude( area );
      newpeak.continuum()->setType( PeakContinuum::Linear );
      newpeak.continuum()->setParameters( chi2abovelinecontstart[index],
                                         chi2abovelinecontcoefs[index],
                                         vector<double>() );
      
      finalpeaks.push_back( newpeak );
    }//while( 1 )
    
    
    //Now need to group peaks together, make them share a ROI and call
    //      std::vector< std::vector<PeakDef> > groupOverlappingPeaks( std::vector<PeakDef> input_peaks );
    const double ncausality = 2.0;  //XXX - arbitrary
    std::vector< std::vector<PeakDef> > groupsofpeaks = causilyDisconnectedPeaks( meas->channel_energies()->front(),
                                                                                 meas->channel_energies()->back(),
                                                                                 ncausality, false, finalpeaks );
    
    finalpeaks.clear();
    for( size_t i = 0; i < groupsofpeaks.size(); ++i )
    {
      double lx = DBL_MAX, ux = -DBL_MAX;
      for( size_t j = 0; j < groupsofpeaks[i].size(); ++j )
      {
        double lowerEnengy, upperEnergy;
        findROIEnergyLimits( lowerEnengy, upperEnergy, groupsofpeaks[i][j], meas, isHPGe );
        lx = std::min( lx, lowerEnengy );
        ux = std::max( ux, upperEnergy );
        if( j )
          groupsofpeaks[i][j].setContinuum( groupsofpeaks[i][0].continuum() );
      }
      
      groupsofpeaks[i][0].continuum()->setRange( lx, ux );
      
      std::shared_ptr<const DetectorPeakResponse> detctorPtr;
      std::vector< std::shared_ptr<const PeakDef> > inputpeaks, resultpeaks;
      
      for( size_t j = 0; j < groupsofpeaks[i].size(); ++j )
      inputpeaks.push_back( std::make_shared<PeakDef>( groupsofpeaks[i][j]) );
      
      
      
      resultpeaks = refitPeaksThatShareROI( meas, detctorPtr, inputpeaks, 0.25 );
      
      for( size_t j = 0; j < resultpeaks.size(); ++j )
      {
        if( (resultpeaks[j]->amplitude()/resultpeaks[j]->amplitudeUncert()) > above_line_chi2_thresh )
          finalpeaks.push_back( *resultpeaks[j] );
      }
    }
    
  }//if( filtered.size() )
  
  
  //now need to filter new peaks that match origpeaks - inefficent
  //  implementation, but whatever for now
  //Note that AutoPeakSearchChi2Fcn already took care of most of the work of not duplicating peaks
  resultpeaks = finalpeaks;
  finalpeaks.clear();
  for( const PeakDef &newpeak : resultpeaks )
  {
    bool matches = false;
    for( const PeakDef &oldpeak : origpeaks )
    {
      matches |= ((fabs(newpeak.mean()-oldpeak.mean())/oldpeak.fwhm()) < 0.1);
    }
    
    if( !matches )
      finalpeaks.push_back( newpeak );
  }//for( size_t i = 0; i < resultpeaks.size(); ++i )
  
  finalpeaks.insert( finalpeaks.end(), origpeaks.begin(), origpeaks.end() );
  std::sort( finalpeaks.begin(), finalpeaks.end(), &PeakDef::lessThanByMean );
  
  
  return finalpeaks;
}//vector<PeakDef> search_for_peaks(...)
}//namespace ExperimentalPeakSearch
        

const char *calculateContinuum( float *spectrum, int ssize,
                               int numberIterations,
                               int direction, int filterOrder,
                               bool smoothing,int smoothWindow,
                               bool compton )
{
  
  //Original Author: Miroslav Morhac   27/05/99
  //Adapted from TSpectrum  class (specifically the Background function)
  //  by wcjohnson 20120216, and slowly being converted to multithreaded code...
  //    basic idea for multithreading:
  //      -first declare all variables only within scope needed (originally all
  //       variable decalared at beinging of this funtion) - Mostly done I think
  //      -break loops up where possible to have each thread work on specific
  //       portion of the array
  //      -Will probably want to use a thread_pool to minimize thread overhead
  //      -Can probably just obptimize the kBackOrder4 portion of the code,
  //       since this seems to work best anyway
  //
  //  See: http://root.cern.ch/root/html/TSpectrum.html for documentation
  //  ROOT code is licenced under the LGPL, see http://root.cern.ch/root/License.html
  
  // April 2012 Notes:
  // Matthew attempted to optimize/parallelize/thread'ize this function, but
  // couldn't implement anything that looked like it'd be faster than the serial
  // version below.  It's not embarassing parallel, and attempts to
  // parallelize or optimize the parts that are had too much overhead and
  // ended up slowing things down.  You can continue where Matthew left
  // off by looking at the following places:
  // (1) The users/branches/calculateContinuumMt branch, which contains
  //     code that tried to parallelize the "for" loops.  This slowed things
  //     down a lot.
  // (2) Revision -r 6768 in the svn trunk, which contains an attempt to
  //     add up values in the working_space only once, instead of up to 7
  //     times. This was slower by about 4x, possible due to sorting the
  //     fenceposts using insertion sort.  Perhaps a sorting network of
  //     size 14 or 16 would be faster, but I couldn't find an implementation
  // There is a "test" called calculateContinuumTest that might be useful
  // for timing tests and verification of possible data corruption.
  
  if (ssize <= 0)
    return "Wrong Parameters";
  if (numberIterations < 1)
    return "Width of Clipping Window Must Be Positive";
  if (ssize < 2 * numberIterations + 1)
    return "Too Large Clipping Window";
  if (smoothing == true && smoothWindow != kBackSmoothing3 && smoothWindow != kBackSmoothing5 && smoothWindow != kBackSmoothing7 && smoothWindow != kBackSmoothing9 && smoothWindow != kBackSmoothing11 && smoothWindow != kBackSmoothing13 && smoothWindow != kBackSmoothing15)
    return "Incorrect width of smoothing window";
  
  
  float *working_space = new float[2 * ssize];
  std::unique_ptr<float []> working_space_scoper( working_space );
  
  for( int i = 0; i < ssize; i++ )
  {
    working_space[i] = spectrum[i];
    working_space[i + ssize] = spectrum[i];
  }//for (i = 0; i < ssize; i++)
  
  const int bw=(smoothWindow-1)/2;
  
  int n_iters = (direction == kBackIncreasingWindow) ? 1 : numberIterations;
  
  if (filterOrder == kBackOrder2)
  {
    do
    {
      for ( int j = n_iters; j < ssize - n_iters; j++)
      {
        if (smoothing == false)
        {
          float a = working_space[ssize + j];
          float b = (working_space[ssize + j - n_iters] + working_space[ssize + j + n_iters]) / 2.0;
          if (b < a)
            a = b;
          working_space[j] = a;
        }else //if (smoothing == true)
        {
          float a = working_space[ssize + j];
          float av = 0;
          float men = 0;
          for ( int w = j - bw; w <= j + bw; w++)
          {
            if ( w >= 0 && w < ssize)
            {
              av += working_space[ssize + w];
              men +=1;
            }
          }
          av = av / men;
          float b = 0;
          men = 0;
          for ( int w = j - n_iters - bw; w <= j - n_iters + bw; w++){
            if ( w >= 0 && w < ssize){
              b += working_space[ssize + w];
              men +=1;
            }
          }
          b = b / men;
          float c = 0;
          men = 0;
          for ( int w = j + n_iters - bw; w <= j + n_iters + bw; w++)
          {
            if ( w >= 0 && w < ssize)
            {
              c += working_space[ssize + w];
              men +=1;
            }
          }
          c = c / men;
          b = (b + c) / 2;
          if (b < a)
            av = b;
          working_space[j]=av;
        }//if (smoothing == false) / else
      }//for (j = n_iters; j < ssize - n_iters; j++)
      
      for (int j = n_iters; j < ssize - n_iters; j++)
      working_space[ssize + j] = working_space[j];
      
      if (direction == kBackIncreasingWindow)
        n_iters+=1;
      else if(direction == kBackDecreasingWindow)
        n_iters-=1;
    }while( (direction == kBackIncreasingWindow && n_iters <= numberIterations)
           || (direction == kBackDecreasingWindow && n_iters >= 1) );
  }//if (filterOrder == kBackOrder2)
  
  else if (filterOrder == kBackOrder4)
  {
    cerr << "kBackOrder4" << endl;
    do{
      for (int j = n_iters; j < ssize - n_iters; j++) {
        if (smoothing == false){
          float a = working_space[ssize + j];
          float b = (working_space[ssize + j - n_iters] + working_space[ssize + j + n_iters]) / 2.0;
          float c = 0;
          float ai = n_iters / 2;
          c -= working_space[ssize + j - (int) (2 * ai)] / 6;
          c += 4 * working_space[ssize + j - (int) ai] / 6;
          c += 4 * working_space[ssize + j + (int) ai] / 6;
          c -= working_space[ssize + j + (int) (2 * ai)] / 6;
          if (b < c)
            b = c;
          if (b < a)
            a = b;
          working_space[j] = a;
        }
        else if (smoothing == true){
          float a = working_space[ssize + j];
          float ai = n_iters / 2;
          float av = 0;
          float men = 0;
          
          for ( int w = j - bw; w <= j + bw; w++){
            if ( w >= 0 && w < ssize){
              av += working_space[ssize + w];
              men +=1;
            }
          }
          av = av / men;
          
          
          float b = 0;
          men = 0;
          for ( int w = j - n_iters - bw; w <= j - n_iters + bw; w++){
            if ( w >= 0 && w < ssize){
              b += working_space[ssize + w];
              men +=1;
            }
          }
          b = b / men;
          
          float c = 0;
          men = 0;
          for ( int w = j + n_iters - bw; w <= j + n_iters + bw; w++){
            if ( w >= 0 && w < ssize){
              c += working_space[ssize + w];
              men +=1;
            }
          }
          c = c / men;
          b = (b + c) / 2;
          
          
          float b4 = 0;
          men = 0;
          for ( int w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              b4 += working_space[ssize + w];
              men +=1;
            }
          }
          b4 = b4 / men;
          
          
          float c4 = 0;
          men = 0;
          for ( int w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              c4 += working_space[ssize + w];
              men +=1;
            }
          }
          c4 = c4 / men;
          
          
          float d4 = 0;
          men = 0;
          for ( int w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              d4 += working_space[ssize + w];
              men +=1;
            }
          }
          d4 = d4 / men;
          
          float e4 = 0;
          men = 0;
          for ( int w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              e4 += working_space[ssize + w];
              men +=1;
            }
          }
          e4 = e4 / men;
          b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
          if (b < b4)
            b = b4;
          if (b < a)
            av = b;
          working_space[j]=av;
        }
      }
      for (int j = n_iters; j < ssize - n_iters; j++)
      working_space[ssize + j] = working_space[j];
      if (direction == kBackIncreasingWindow)
        n_iters+=1;
      else if(direction == kBackDecreasingWindow)
        n_iters-=1;
    }while((direction == kBackIncreasingWindow && n_iters <= numberIterations) || (direction == kBackDecreasingWindow && n_iters >= 1));
  }
  
  else if (filterOrder == kBackOrder6) {
    do{
      for (int j = n_iters; j < ssize - n_iters; j++) {
        
        if (smoothing == false){
          float a = working_space[ssize + j];
          float b = (working_space[ssize + j - n_iters] + working_space[ssize + j + n_iters]) / 2.0;
          float c = 0;
          float ai = n_iters / 2;
          c -= working_space[ssize + j - (int) (2 * ai)] / 6;
          c += 4 * working_space[ssize + j - (int) ai] / 6;
          c += 4 * working_space[ssize + j + (int) ai] / 6;
          c -= working_space[ssize + j + (int) (2 * ai)] / 6;
          float d = 0;
          ai = n_iters / 3;
          d += working_space[ssize + j - (int) (3 * ai)] / 20;
          d -= 6 * working_space[ssize + j - (int) (2 * ai)] / 20;
          d += 15 * working_space[ssize + j - (int) ai] / 20;
          d += 15 * working_space[ssize + j + (int) ai] / 20;
          d -= 6 * working_space[ssize + j + (int) (2 * ai)] / 20;
          d += working_space[ssize + j + (int) (3 * ai)] / 20;
          if (b < d)
            b = d;
          if (b < c)
            b = c;
          if (b < a)
            a = b;
          working_space[j] = a;
        }
        
        else if (smoothing == true){
          float a = working_space[ssize + j];
          float av = 0;
          float men = 0;
          for ( int w = j - bw; w <= j + bw; w++){
            if ( w >= 0 && w < ssize){
              av += working_space[ssize + w];
              men +=1;
            }
          }
          av = av / men;
          float b = 0;
          men = 0;
          for ( int w = j - n_iters - bw; w <= j - n_iters + bw; w++){
            if ( w >= 0 && w < ssize){
              b += working_space[ssize + w];
              men +=1;
            }
          }
          b = b / men;
          float c = 0;
          men = 0;
          for ( int w = j + n_iters - bw; w <= j + n_iters + bw; w++){
            if ( w >= 0 && w < ssize){
              c += working_space[ssize + w];
              men +=1;
            }
          }
          c = c / men;
          b = (b + c) / 2;
          float ai = n_iters / 2;
          float b4 = 0;
          men = 0;
          for ( int w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              b4 += working_space[ssize + w];
              men +=1;
            }
          }
          b4 = b4 / men;
          float c4 = 0;
          men = 0;
          for ( int w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              c4 += working_space[ssize + w];
              men +=1;
            }
          }
          c4 = c4 / men;
          float d4 = 0;
          men = 0;
          for ( int w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              d4 += working_space[ssize + w];
              men +=1;
            }
          }
          d4 = d4 / men;
          float e4 = 0;
          men = 0;
          for ( int w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              e4 += working_space[ssize + w];
              men +=1;
            }
          }
          e4 = e4 / men;
          b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
          ai = n_iters / 3;
          float b6 = 0;
          men = 0;
          for ( int w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              b6 += working_space[ssize + w];
              men +=1;
            }
          }
          b6 = b6 / men;
          float c6 = 0;
          men = 0;
          for ( int w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              c6 += working_space[ssize + w];
              men +=1;
            }
          }
          c6 = c6 / men;
          float d6 = 0;
          men = 0;
          for ( int w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              d6 += working_space[ssize + w];
              men +=1;
            }
          }
          d6 = d6 / men;
          float e6 = 0;
          men = 0;
          for ( int w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              e6 += working_space[ssize + w];
              men +=1;
            }
          }
          e6 = e6 / men;
          float f6 = 0;
          men = 0;
          for ( int w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              f6 += working_space[ssize + w];
              men +=1;
            }
          }
          f6 = f6 / men;
          float g6 = 0;
          men = 0;
          for ( int w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              g6 += working_space[ssize + w];
              men +=1;
            }
          }
          g6 = g6 / men;
          b6 = (b6 - 6 * c6 + 15 * d6 + 15 * e6 - 6 * f6 + g6) / 20;
          if (b < b6)
            b = b6;
          if (b < b4)
            b = b4;
          if (b < a)
            av = b;
          working_space[j]=av;
        }
      }
      for (int j = n_iters; j < ssize - n_iters; j++)
      working_space[ssize + j] = working_space[j];
      if (direction == kBackIncreasingWindow)
        n_iters+=1;
      else if(direction == kBackDecreasingWindow)
        n_iters-=1;
    }while((direction == kBackIncreasingWindow && n_iters <= numberIterations) || (direction == kBackDecreasingWindow && n_iters >= 1));
  }
  
  else if (filterOrder == kBackOrder8) {
    do{
      for (int j = n_iters; j < ssize - n_iters; j++) {
        if (smoothing == false){
          float a = working_space[ssize + j];
          float b = (working_space[ssize + j - n_iters] + working_space[ssize + j + n_iters]) / 2.0;
          float c = 0;
          float ai = n_iters / 2;
          c -= working_space[ssize + j - (int) (2 * ai)] / 6;
          c += 4 * working_space[ssize + j - (int) ai] / 6;
          c += 4 * working_space[ssize + j + (int) ai] / 6;
          c -= working_space[ssize + j + (int) (2 * ai)] / 6;
          float d = 0;
          ai = n_iters / 3;
          d += working_space[ssize + j - (int) (3 * ai)] / 20;
          d -= 6 * working_space[ssize + j - (int) (2 * ai)] / 20;
          d += 15 * working_space[ssize + j - (int) ai] / 20;
          d += 15 * working_space[ssize + j + (int) ai] / 20;
          d -= 6 * working_space[ssize + j + (int) (2 * ai)] / 20;
          d += working_space[ssize + j + (int) (3 * ai)] / 20;
          float e = 0;
          ai = n_iters / 4;
          e -= working_space[ssize + j - (int) (4 * ai)] / 70;
          e += 8 * working_space[ssize + j - (int) (3 * ai)] / 70;
          e -= 28 * working_space[ssize + j - (int) (2 * ai)] / 70;
          e += 56 * working_space[ssize + j - (int) ai] / 70;
          e += 56 * working_space[ssize + j + (int) ai] / 70;
          e -= 28 * working_space[ssize + j + (int) (2 * ai)] / 70;
          e += 8 * working_space[ssize + j + (int) (3 * ai)] / 70;
          e -= working_space[ssize + j + (int) (4 * ai)] / 70;
          if (b < e)
            b = e;
          if (b < d)
            b = d;
          if (b < c)
            b = c;
          if (b < a)
            a = b;
          working_space[j] = a;
        }
        
        else if (smoothing == true)
        {
          float a = working_space[ssize + j];
          float av = 0;
          float men = 0;
          
          for ( int w = j - bw; w <= j + bw; w++){
            if ( w >= 0 && w < ssize){
              av += working_space[ssize + w];
              men +=1;
            }
          }
          av = av / men;
          float b = 0;
          men = 0;
          for ( int w = j - n_iters - bw; w <= j - n_iters + bw; w++){
            if ( w >= 0 && w < ssize){
              b += working_space[ssize + w];
              men +=1;
            }
          }
          b = b / men;
          float c = 0;
          men = 0;
          for ( int w = j + n_iters - bw; w <= j + n_iters + bw; w++){
            if ( w >= 0 && w < ssize){
              c += working_space[ssize + w];
              men +=1;
            }
          }
          c = c / men;
          b = (b + c) / 2;
          float ai = n_iters / 2;
          float b4 = 0;
          men = 0;
          for ( int w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              b4 += working_space[ssize + w];
              men +=1;
            }
          }
          b4 = b4 / men;
          float c4 = 0;
          men = 0;
          for ( int w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              c4 += working_space[ssize + w];
              men +=1;
            }
          }
          c4 = c4 / men;
          float d4 = 0;
          men = 0;
          for ( int w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              d4 += working_space[ssize + w];
              men +=1;
            }
          }
          d4 = d4 / men;
          float e4 = 0;
          men = 0;
          for ( int w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              e4 += working_space[ssize + w];
              men +=1;
            }
          }
          e4 = e4 / men;
          b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
          ai = n_iters / 3;
          float b6 = 0;
          men = 0;
          for ( int w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              b6 += working_space[ssize + w];
              men +=1;
            }
          }
          b6 = b6 / men;
          float c6 = 0;
          men = 0;
          for ( int w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              c6 += working_space[ssize + w];
              men +=1;
            }
          }
          c6 = c6 / men;
          float d6 = 0;
          men = 0;
          for ( int w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              d6 += working_space[ssize + w];
              men +=1;
            }
          }
          d6 = d6 / men;
          float e6 = 0;
          men = 0;
          for ( int w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              e6 += working_space[ssize + w];
              men +=1;
            }
          }
          e6 = e6 / men;
          float f6 = 0;
          men = 0;
          for ( int w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              f6 += working_space[ssize + w];
              men +=1;
            }
          }
          f6 = f6 / men;
          float g6 = 0;
          men = 0;
          for ( int w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              g6 += working_space[ssize + w];
              men +=1;
            }
          }
          g6 = g6 / men;
          b6 = (b6 - 6 * c6 + 15 * d6 + 15 * e6 - 6 * f6 + g6) / 20;
          ai = n_iters / 4;
          float b8 = 0;
          men = 0;
          for ( int w = j - (int)(4 * ai) - bw; w <= j - (int)(4 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              b8 += working_space[ssize + w];
              men +=1;
            }
          }
          b8 = b8 / men;
          float c8 = 0;
          men = 0;
          for ( int w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              c8 += working_space[ssize + w];
              men +=1;
            }
          }
          c8 = c8 / men;
          float d8 = 0;
          men = 0;
          for ( int w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              d8 += working_space[ssize + w];
              men +=1;
            }
          }
          d8 = d8 / men;
          float e8 = 0;
          men = 0;
          for ( int w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              e8 += working_space[ssize + w];
              men +=1;
            }
          }
          e8 = e8 / men;
          float f8 = 0;
          men = 0;
          for ( int w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
            if (w >= 0 && w < ssize){
              f8 += working_space[ssize + w];
              men +=1;
            }
          }
          f8 = f8 / men;
          float g8 = 0;
          men = 0;
          for ( int w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              g8 += working_space[ssize + w];
              men +=1;
            }
          }
          g8 = g8 / men;
          float h8 = 0;
          men = 0;
          for ( int w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              h8 += working_space[ssize + w];
              men +=1;
            }
          }
          h8 = h8 / men;
          float i8 = 0;
          men = 0;
          for ( int w = j + (int)(4 * ai) - bw; w <= j + (int)(4 * ai) + bw; w++){
            if (w >= 0 && w < ssize){
              i8 += working_space[ssize + w];
              men +=1;
            }
          }
          i8 = i8 / men;
          b8 = ( -b8 + 8 * c8 - 28 * d8 + 56 * e8 - 56 * f8 - 28 * g8 + 8 * h8 - i8)/70;
          if (b < b8)
            b = b8;
          if (b < b6)
            b = b6;
          if (b < b4)
            b = b4;
          if (b < a)
            av = b;
          working_space[j]=av;
        }
      }
      for (int j = n_iters; j < ssize - n_iters; j++)
      working_space[ssize + j] = working_space[j];
      if (direction == kBackIncreasingWindow)
        n_iters += 1;
      else if(direction == kBackDecreasingWindow)
        n_iters -= 1;
    }while((direction == kBackIncreasingWindow && n_iters <= numberIterations) || (direction == kBackDecreasingWindow && n_iters >= 1));
  }
  
  if (compton == true) {
    int b2 = 0;
    for (int i = 0; i < ssize; i++){
      int b1 = b2;
      float a = working_space[i], b = spectrum[i];
      
      //           j = i;
      
      if (fabs(a - b) >= 1) {
        b1 = i - 1;
        if (b1 < 0)
          b1 = 0;
        float yb1 = working_space[b1];
        float c = 0.0;
        int priz = 0;
        for (b2 = b1 + 1; priz == 0 && b2 < ssize; b2++){
          a = working_space[b2], b = spectrum[b2];
          c = c + b - yb1;
          if (fabs(a - b) < 1) {
            priz = 1;
            //                    yb2 = b;
          }
        }
        if (b2 == ssize)
          b2 -= 1;
        
        float yb2 = working_space[b2];
        
        if (yb1 <= yb2)
        {
          c = 0.0;
          for (int j = b1; j <= b2; j++){
            b = spectrum[j];
            c = c + b - yb1;
          }
          if (c > 1){
            c = (yb2 - yb1) / c;
            float d = 0.0;
            for (int j = b1; j <= b2 && j < ssize; j++){
              b = spectrum[j];
              d = d + b - yb1;
              a = c * d + yb1;
              working_space[ssize + j] = a;
            }
          }
        }else
        {
          c = 0.0;
          for (int j = b2; j >= b1; j--){
            b = spectrum[j];
            c = c + b - yb2;
          }
          if (c > 1){
            c = (yb1 - yb2) / c;
            float d = 0.0;
            for (int j = b2;j >= b1 && j >= 0; j--){
              b = spectrum[j];
              d = d + b - yb2;
              a = c * d + yb2;
              working_space[ssize + j] = a;
            }
          }
        }
        i=b2;
      }
    }
  }
  
  for (int j = 0; j < ssize; j++)
  spectrum[j] = working_space[ssize + j];
  //     delete[]working_space;
  
  /*
   boost::posix_time::ptime endTime = boost::posix_time::microsec_clock::local_time();
   boost::posix_time::time_duration durationEnd( endTime.time_of_day() );
   
   //   std::cout << "start duration is " << durationStart.total_milliseconds() << std::endl;
   //   std::cout << "end duration is " << durationEnd.total_milliseconds() << std::endl;
   
   boost::posix_time::time_duration duration = durationEnd - durationStart;
   double seconds = duration.total_milliseconds()/1000.0;
   std::cout << "****duration is " << seconds << " seconds" << std::endl;
   string command;
   stringstream out;
   out << "banner " << seconds;
   command = out.str();
   system( command.c_str() );
   */
  
  return 0;
}//const char *calculateContinuum(...)







vector<float> findPeaksByRelaxation( float *source,float *destVector, int ssize,
                                    float sigma, double threshold,
                                    bool backgroundRemove,int deconIterations,
                                    bool markov, int averWindow )
{
  //Author: Miroslav Morhac   27/05/99
  //Adapted from TSpectrum  class (specifically Background)
  //  by wcjohnson 20120216.
  //  See: http://root.cern.ch/root/html/TSpectrum.html for documentation
  
#define PEAK_WINDOW 1024
  const int fMaxPeaks = 100;
  vector<float> fPositionX;
  
  int i, j, numberIterations = (int)(7 * sigma + 0.5);
  double a, b, c;
  int k, lindex, posit, imin, imax, jmin, jmax, lh_gold, priz;
  double lda, ldb, ldc, area, maximum, maximum_decon;
  int xmin, xmax, l, peak_index = 0, size_ext = ssize + 2 * numberIterations, shift = numberIterations, bw = 2;
  
  double maxch;
  double nom, nip, nim, sp, sm, plocha = 0;
  double m0low=0,m1low=0,m2low=0,l0low=0,l1low=0,detlow,av,men;
  if (sigma < 1) {
    cerr << "findPeaksByRelaxation(...)\n\t"
    << "Invalid sigma, must be greater than or equal to 1" << endl;
    return fPositionX;
  }
  
  if(threshold<=0 || threshold>=100){
    cerr << "findPeaksByRelaxation(...)\n\t"
    << "Invalid threshold, must be positive and less than 100"
    << endl;
    return fPositionX;
  }
  
  j = (int) (5.0 * sigma + 0.5);
  if (j >= PEAK_WINDOW / 2) {
    cerr << "findPeaksByRelaxation(...)\n\tToo large sigma" << endl;
    return fPositionX;
  }
  
  if (markov == true) {
    if (averWindow <= 0) {
      cerr << "findPeaksByRelaxation(...)\n\t"
      << "Averanging window must be positive" << endl;
      return fPositionX;
    }
  }
  
  if(backgroundRemove == true){
    if(ssize < 2 * numberIterations + 1){
      cerr << "findPeaksByRelaxation(...)\n\t"
      << "Too large clipping window" << endl;
      return fPositionX;
    }
  }
  
  fPositionX.resize( fMaxPeaks, 0.0 );
  
  k = int(2 * sigma+0.5);
  if(k >= 2){
    for(i = 0;i < k;i++){
      a = i,b = source[i];
      m0low += 1,m1low += a,m2low += a * a,l0low += b,l1low += a * b;
    }
    detlow = m0low * m2low - m1low * m1low;
    if(detlow != 0)
      l1low = (-l0low * m1low + l1low * m0low) / detlow;
    
    else
      l1low = 0;
    if(l1low > 0)
      l1low=0;
  }
  
  else{
    l1low = 0;
  }
  
  i = (int)(7 * sigma + 0.5);
  i = 2 * i;
  double *working_space = new double [7 * (ssize + i)];
  for (j=0;j<7 * (ssize + i);j++) working_space[j] = 0;
  for(i = 0; i < size_ext; i++){
    if(i < shift){
      a = i - shift;
      working_space[i + size_ext] = source[0] + l1low * a;
      if(working_space[i + size_ext] < 0)
        working_space[i + size_ext]=0;
    }
    
    else if(i >= ssize + shift){
      a = i - (ssize - 1 + shift);
      working_space[i + size_ext] = source[ssize - 1];
      if(working_space[i + size_ext] < 0)
        working_space[i + size_ext]=0;
    }
    
    else
      working_space[i + size_ext] = source[i - shift];
  }
  
  if(backgroundRemove == true){
    for(i = 1; i <= numberIterations; i++){
      for(j = i; j < size_ext - i; j++){
        if(markov == false){
          a = working_space[size_ext + j];
          b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
          if(b < a)
            a = b;
          
          working_space[j]=a;
        }
        
        else{
          a = working_space[size_ext + j];
          av = 0;
          men = 0;
          for ( int w = j - bw; w <= j + bw; w++){
            if ( w >= 0 && w < size_ext){
              av += working_space[size_ext + w];
              men +=1;
            }
          }
          av = av / men;
          b = 0;
          men = 0;
          for ( int w = j - i - bw; w <= j - i + bw; w++){
            if ( w >= 0 && w < size_ext){
              b += working_space[size_ext + w];
              men +=1;
            }
          }
          b = b / men;
          c = 0;
          men = 0;
          for ( int w = j + i - bw; w <= j + i + bw; w++){
            if ( w >= 0 && w < size_ext){
              c += working_space[size_ext + w];
              men +=1;
            }
          }
          c = c / men;
          b = (b + c) / 2;
          if (b < a)
            av = b;
          working_space[j]=av;
        }
      }
      for(j = i; j < size_ext - i; j++)
      working_space[size_ext + j] = working_space[j];
    }
    for(j = 0;j < size_ext; j++){
      if(j < shift){
        a = j - shift;
        b = source[0] + l1low * a;
        if (b < 0) b = 0;
        working_space[size_ext + j] = b - working_space[size_ext + j];
      }
      
      else if(j >= ssize + shift){
        a = j - (ssize - 1 + shift);
        b = source[ssize - 1];
        if (b < 0) b = 0;
        working_space[size_ext + j] = b - working_space[size_ext + j];
      }
      
      else{
        working_space[size_ext + j] = source[j - shift] - working_space[size_ext + j];
      }
    }
    for(j = 0;j < size_ext; j++){
      if(working_space[size_ext + j] < 0) working_space[size_ext + j] = 0;
    }
  }
  
  for(i = 0; i < size_ext; i++){
    working_space[i + 6*size_ext] = working_space[i + size_ext];
  }
  
  if(markov == true){
    for(j = 0; j < size_ext; j++)
    working_space[2 * size_ext + j] = working_space[size_ext + j];
    xmin = 0,xmax = size_ext - 1;
    for(i = 0, maxch = 0; i < size_ext; i++){
      working_space[i] = 0;
      if(maxch < working_space[2 * size_ext + i])
        maxch = working_space[2 * size_ext + i];
      plocha += working_space[2 * size_ext + i];
    }
    if(maxch == 0) {
      delete [] working_space;
      fPositionX.clear();
      return fPositionX;;
    }
    
    nom = 1;
    working_space[xmin] = 1;
    for(i = xmin; i < xmax; i++){
      nip = working_space[2 * size_ext + i] / maxch;
      nim = working_space[2 * size_ext + i + 1] / maxch;
      sp = 0,sm = 0;
      for(l = 1; l <= averWindow; l++){
        if((i + l) > xmax)
          a = working_space[2 * size_ext + xmax] / maxch;
        
        else
          a = working_space[2 * size_ext + i + l] / maxch;
        
        b = a - nip;
        if(a + nip <= 0)
          a=1;
        
        else
          a = sqrt(a + nip);
        
        b = b / a;
        b = exp(b);
        sp = sp + b;
        if((i - l + 1) < xmin)
          a = working_space[2 * size_ext + xmin] / maxch;
        
        else
          a = working_space[2 * size_ext + i - l + 1] / maxch;
        
        b = a - nim;
        if(a + nim <= 0)
          a = 1;
        
        else
          a = sqrt(a + nim);
        
        b = b / a;
        b = exp(b);
        sm = sm + b;
      }
      a = sp / sm;
      a = working_space[i + 1] = working_space[i] * a;
      nom = nom + a;
    }
    for(i = xmin; i <= xmax; i++){
      working_space[i] = working_space[i] / nom;
    }
    for(j = 0; j < size_ext; j++)
    working_space[size_ext + j] = working_space[j] * plocha;
    for(j = 0; j < size_ext; j++){
      working_space[2 * size_ext + j] = working_space[size_ext + j];
    }
    if(backgroundRemove == true){
      for(i = 1; i <= numberIterations; i++){
        for(j = i; j < size_ext - i; j++){
          a = working_space[size_ext + j];
          b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
          if(b < a)
            a = b;
          working_space[j] = a;
        }
        for(j = i; j < size_ext - i; j++)
        working_space[size_ext + j] = working_space[j];
      }
      for(j = 0; j < size_ext; j++){
        working_space[size_ext + j] = working_space[2 * size_ext + j] - working_space[size_ext + j];
      }
    }
  }
  //deconvolution starts
  area = 0;
  lh_gold = -1;
  posit = 0;
  maximum = 0;
  //generate response vector
  for(i = 0; i < size_ext; i++){
    lda = (double)i - 3 * sigma;
    lda = lda * lda / (2 * sigma * sigma);
    j = (int)(1000 * exp(-lda));
    lda = j;
    if(lda != 0)
      lh_gold = i + 1;
    
    working_space[i] = lda;
    area = area + lda;
    if(lda > maximum){
      maximum = lda;
      posit = i;
    }
  }
  //read source vector
  for(i = 0; i < size_ext; i++)
  working_space[2 * size_ext + i] = fabs(working_space[size_ext + i]);
  //create matrix at*a(vector b)
  i = lh_gold - 1;
  if(i > size_ext)
    i = size_ext;
  
  imin = -i,imax = i;
  for(i = imin; i <= imax; i++){
    lda = 0;
    jmin = 0;
    if(i < 0)
      jmin = -i;
    jmax = lh_gold - 1 - i;
    if(jmax > (lh_gold - 1))
      jmax = lh_gold - 1;
    
    for(j = jmin;j <= jmax; j++){
      ldb = working_space[j];
      ldc = working_space[i + j];
      lda = lda + ldb * ldc;
    }
    working_space[size_ext + i - imin] = lda;
  }
  //create vector p
  i = lh_gold - 1;
  imin = -i,imax = size_ext + i - 1;
  for(i = imin; i <= imax; i++){
    lda = 0;
    for(j = 0; j <= (lh_gold - 1); j++){
      ldb = working_space[j];
      k = i + j;
      if(k >= 0 && k < size_ext){
        ldc = working_space[2 * size_ext + k];
        lda = lda + ldb * ldc;
      }
      
    }
    working_space[4 * size_ext + i - imin] = lda;
  }
  //move vector p
  for(i = imin; i <= imax; i++)
  working_space[2 * size_ext + i - imin] = working_space[4 * size_ext + i - imin];
  //initialization of resulting vector
  for(i = 0; i < size_ext; i++)
  working_space[i] = 1;
  //START OF ITERATIONS
  for(lindex = 0; lindex < deconIterations; lindex++){
    for(i = 0; i < size_ext; i++){
      if(fabs(working_space[2 * size_ext + i]) > 0.00001 && fabs(working_space[i]) > 0.00001){
        lda=0;
        jmin = lh_gold - 1;
        if(jmin > i)
          jmin = i;
        
        jmin = -jmin;
        jmax = lh_gold - 1;
        if(jmax > (size_ext - 1 - i))
          jmax=size_ext-1-i;
        
        for(j = jmin; j <= jmax; j++){
          ldb = working_space[j + lh_gold - 1 + size_ext];
          ldc = working_space[i + j];
          lda = lda + ldb * ldc;
        }
        ldb = working_space[2 * size_ext + i];
        if(lda != 0)
          lda = ldb / lda;
        
        else
          lda = 0;
        
        ldb = working_space[i];
        lda = lda * ldb;
        working_space[3 * size_ext + i] = lda;
      }
    }
    for(i = 0; i < size_ext; i++){
      working_space[i] = working_space[3 * size_ext + i];
    }
  }
  //shift resulting spectrum
  for(i=0;i<size_ext;i++){
    lda = working_space[i];
    j = i + posit;
    j = j % size_ext;
    working_space[size_ext + j] = lda;
  }
  //write back resulting spectrum
  maximum = 0, maximum_decon = 0;
  j = lh_gold - 1;
  for(i = 0; i < size_ext - j; i++){
    if(i >= shift && i < ssize + shift){
      working_space[i] = area * working_space[size_ext + i + j];
      if(maximum_decon < working_space[i])
        maximum_decon = working_space[i];
      if(maximum < working_space[6 * size_ext + i])
        maximum = working_space[6 * size_ext + i];
    }
    
    else
      working_space[i] = 0;
  }
  lda=1;
  if(lda>threshold)
    lda=threshold;
  lda=lda/100;
  
  //searching for peaks in deconvolved spectrum
  for(i = 1; i < size_ext - 1; i++){
    if(working_space[i] > working_space[i - 1] && working_space[i] > working_space[i + 1]){
      if(i >= shift && i < ssize + shift){
        if(working_space[i] > lda*maximum_decon && working_space[6 * size_ext + i] > threshold * maximum / 100.0){
          for(j = i - 1, a = 0, b = 0; j <= i + 1; j++){
            a += (double)(j - shift) * working_space[j];
            b += working_space[j];
          }
          a = a / b;
          if(a < 0)
            a = 0;
          
          if(a >= ssize)
            a = ssize - 1;
          if(peak_index == 0){
            fPositionX[0] = a;
            peak_index = 1;
          }
          
          else{
            for(j = 0, priz = 0; j < peak_index && priz == 0; j++){
              if(working_space[6 * size_ext + shift + (int)a] > working_space[6 * size_ext + shift + (int)fPositionX[j]])
                priz = 1;
            }
            if(priz == 0){
              if(j < fMaxPeaks){
                fPositionX[j] = a;
              }
            }
            
            else{
              for(k = peak_index; k >= j; k--){
                if(k < fMaxPeaks){
                  fPositionX[k] = fPositionX[k - 1];
                }
              }
              fPositionX[j - 1] = a;
            }
            if(peak_index < fMaxPeaks)
              peak_index += 1;
          }
        }
      }
    }
  }
  
  for(i = 0; i < ssize; i++) destVector[i] = working_space[i + shift];
  delete [] working_space;
  
  if(peak_index == fMaxPeaks)
    cerr << "findPeaksByRelaxation(...)\n\tPeak buffer full" << endl;
  
  fPositionX.resize( peak_index );
  
  return fPositionX;
}//vector<float> findPeaksByRelaxation(...)








        
