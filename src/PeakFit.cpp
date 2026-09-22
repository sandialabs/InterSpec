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
#include "Eigen/Dense"



#include "SpecUtils/SpecFile.h"
#include "SpecUtils/SpecUtilsAsync.h"
//#include "SpecUtils/DateTime.h" //on for debug SpecUtils::get_cpu_time() / SpecUtils::get_wall_time();

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitLM.h"
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/DetectionLimitCalc.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/DetectorPeakResponse.h"



#include "InterSpec/PeakFit_imp.hpp"


using namespace std;

using SpecUtils::Measurement;


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

namespace
{
  /** Cooperative cancellation check used throughout the peak-search call graph.
   `nullptr` token == never canceled.  Loaded with relaxed ordering since the
   only thing the reader does on success is bail out of a loop.
   */
  inline bool is_search_canceled( const std::shared_ptr<const std::atomic<bool>> &f )
  {
    return f && f->load( std::memory_order_relaxed );
  }
}//namespace

bool largerByAmplitude( const std::shared_ptr<const PeakDef> &lhs, const std::shared_ptr<const PeakDef> &rhs )
{
  return lhs->amplitude() > rhs->amplitude();
}


void do_peak_automated_searchfit( const double x,
                                  const std::shared_ptr<const Measurement> &meas,
                                  const std::shared_ptr<const DetectorPeakResponse> &drf,
                                  const PeakShrdVec &inpeaks,
                                  std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                                  std::pair< PeakShrdVec, PeakShrdVec > &answer,
                                  const std::shared_ptr<const std::atomic<bool>> &cancel_flag )
{
  if( is_search_canceled( cancel_flag ) )
  {
    answer.first.clear();
    answer.second.clear();
    return;
  }

  try
  {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
    DebugLog(cout) << "Will try fitting peak clicked on at " << x << "\n";
#endif
    answer = searchForPeakFromUser( x, -1.0, meas, inpeaks, drf, nullptr, fitPrefs, cancel_flag );
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
                          std::vector<std::shared_ptr<const PeakDef> > input,
                          std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                          const std::shared_ptr<const std::atomic<bool>> &cancel_flag )
{
  if( !meas )
    return input;

  if( is_search_canceled( cancel_flag ) )
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
    if( is_search_canceled( cancel_flag ) )
    {
      // Return whatever we've kept so far, plus the remaining unprocessed peaks
      // un-evaluated.  Cancellation should not destroy data.
      for( size_t j = i; j < npeaks; ++j )
        answer.push_back( input[j] );
      return answer;
    }

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
        twoPeaksPlus = searchForPeakFromUser( m + s, -1.0, meas, onepeak, nullptr, nullptr, fitPrefs, cancel_flag );
        twoPeaksMinus = searchForPeakFromUser( m - s, -1.0, meas, onepeak, nullptr, nullptr, fitPrefs, cancel_flag );
        
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
                                       std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                                       const std::shared_ptr<const std::atomic<bool>> &cancel_flag )
{
  typedef std::shared_ptr<PeakDef> PeakPtr;
  typedef std::shared_ptr<const PeakDef> PeakConstPtr;

  assert( fitPrefs );
  const bool isHPGe = (PeakFitUtils::effective_det_type( fitPrefs, meas, nullptr )
                       == PeakFitUtils::CoarseResolutionType::High);

  size_t lower_channel = 0, upper_channel = 0;
  //    PeakFitUtils::find_spectroscopic_extent( meas, lower_channel, upper_channel );
  //    cout << "Start at " << meas->gamma_channel_center( lower_channel ) << " and going through "
  //    << meas->gamma_channel_center( upper_channel ) << endl;
  
  const vector<std::shared_ptr<PeakDef> > initialcandidates
    = secondDerivativePeakCanidatesWithROI( meas, fitPrefs, lower_channel, upper_channel );

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


  // Sanity cap on the outer loop.  Each iteration should remove all the
  //  `candidatesBeingFitFor` entries from `candidates`, so in normal operation
  //  this completes in well under `initialcandidates.size()` iterations.  The
  //  cap is a defense-in-depth guard against an unexpected logic regression
  //  causing the loop to spin without making progress.
  const size_t max_outer_iterations = 2 * initialcandidates.size() + 16;
  size_t outer_iter = 0;

  while( !candidates.empty() )
  {
    if( is_search_canceled( cancel_flag ) )
      return fitpeakvec;

    if( ++outer_iter > max_outer_iterations )
    {
      cerr << "search_for_peaks_multithread: aborting after " << outer_iter
           << " outer iterations (initial candidates: " << initialcandidates.size()
           << ", remaining: " << candidates.size() << ")" << endl;
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "outer-loop iteration cap exceeded" );
#endif
      break;
    }

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
      pair<PeakShrdVec, PeakShrdVec> &res = results[i];

      pool.post(
        [mean, fitPrefs, &res,
         meas_cref = std::as_const(meas),
         drf_cref = std::as_const(drf),
         fitpeakvec_cref = std::as_const(fitpeakvec),
         cancel_flag
        ](){
          do_peak_automated_searchfit( mean, meas_cref, drf_cref, fitpeakvec_cref, fitPrefs, res, cancel_flag );
        }
      );
    }//for( size_t i = 0; i < candidatesBeingFitFor.size(); ++i )


    pool.join();

    if( is_search_canceled( cancel_flag ) )
      return fitpeakvec;
    
    //const double end_cpu_time = SpecUtils::get_cpu_time();
    //const double end_wall_time = SpecUtils::get_wall_time();
    
    //cout << "Cpu  time: " << (end_cpu_time - start_cpu_time) << "s, vs sum " << sum_cpu_time.load() << "s" << endl;
    //cout << "Wall time: " << (end_wall_time - start_wall_time) << "s, vs sum " << sum_wall_time.load() << "s" << endl;
    
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
  
  if( is_search_canceled( cancel_flag ) )
    return fitpeakvec;

  const auto detResolution = PeakFitUtils::coarse_resolution_from_peaks( fitpeakvec );
  if( detResolution == PeakFitUtils::CoarseResolutionType::High )
    fitpeakvec = filter_anomolous_width_peaks_highres( meas, fitpeakvec, fitPrefs, cancel_flag );

  return fitpeakvec;
}//search_for_peaks_multithread(...)
  
  

vector<std::shared_ptr<const PeakDef> > search_for_peaks_singlethread(
                        const std::shared_ptr<const Measurement> meas,
                        const std::shared_ptr<const DetectorPeakResponse> &drf,
                        std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > origpeaks,
                        std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                        const std::shared_ptr<const std::atomic<bool>> &cancel_flag )
{
  typedef std::shared_ptr<PeakDef> PeakPtr;
  typedef std::shared_ptr<const PeakDef> PeakConstPtr;

  assert( fitPrefs );
  const bool isHPGe = (PeakFitUtils::effective_det_type( fitPrefs, meas, nullptr )
                       == PeakFitUtils::CoarseResolutionType::High);

  size_t lower_channel = 0, upper_channel = 0;
  vector<PeakPtr> candidates
   = secondDerivativePeakCanidatesWithROI( meas, fitPrefs, lower_channel, upper_channel );

  const size_t initial_candidate_count = candidates.size();
  const size_t max_outer_iterations = 2 * initial_candidate_count + 16;
  size_t outer_iter = 0;
  
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
    if( is_search_canceled( cancel_flag ) )
      return fitpeakvec;

    if( ++outer_iter > max_outer_iterations )
    {
      cerr << "search_for_peaks_singlethread: aborting after " << outer_iter
           << " outer iterations (initial candidates: " << initial_candidate_count
           << ", remaining: " << candidates.size() << ")" << endl;
#if( PERFORM_DEVELOPER_CHECKS )
      log_developer_error( __func__, "outer-loop iteration cap exceeded" );
#endif
      break;
    }

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
    do_peak_automated_searchfit( p.mean(), meas, drf, fitpeakvec, fitPrefs, results, cancel_flag );
    
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
  
  if( is_search_canceled( cancel_flag ) )
    return fitpeakvec;

  const auto detResolution = PeakFitUtils::coarse_resolution_from_peaks( fitpeakvec );
  if( detResolution == PeakFitUtils::CoarseResolutionType::High )
    fitpeakvec = filter_anomolous_width_peaks_highres( meas, fitpeakvec, fitPrefs, cancel_flag );

  return fitpeakvec;
}//search_for_peaks_singlethread(...)


vector<std::shared_ptr<const PeakDef> > search_for_peaks(
                              const std::shared_ptr<const Measurement> meas,
                              const std::shared_ptr<const DetectorPeakResponse> drf,
                              std::shared_ptr<const deque< std::shared_ptr<const PeakDef> > > origpeaks,
                              const bool singleThreaded,
                              std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                              std::shared_ptr<const std::atomic<bool>> cancel_flag )
{
  vector<std::shared_ptr<const PeakDef> > answer;

  // Library/test callers do not necessarily have application- or DRF-owned
  // preferences.  The lower-level search routines require a non-null object,
  // so establish the documented fallback once at this public entry point.
  // Inferring the detector class here also makes the behavior independent of
  // assertions: Debug and Release now take the same path.
  if( !fitPrefs )
  {
    auto defaults = std::make_shared<PeakFitDetPrefs>();
    defaults->m_det_type = PeakFitUtils::effective_det_type( nullptr, meas, nullptr );
    defaults->m_source = PeakFitDetPrefs::LoadingSource::FromSpectralData;
    fitPrefs = std::move(defaults);
  }

  if( singleThreaded )
    answer = search_for_peaks_singlethread( meas, drf, origpeaks, fitPrefs, cancel_flag );
  else
    answer = search_for_peaks_multithread( meas, drf, origpeaks, fitPrefs, cancel_flag );

  return answer;
}


std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>>
recover_background_peaks_under_foreground(
    const std::vector<std::shared_ptr<const PeakDef>> &foreground_peaks,
    const std::shared_ptr<const SpecUtils::Measurement> &background_spectrum,
    const std::shared_ptr<const std::deque<std::shared_ptr<const PeakDef>>> &background_auto_peaks,
    const std::shared_ptr<const DetectorPeakResponse> &background_drf,
    std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
    std::shared_ptr<const std::atomic<bool>> cancel_flag )
{
  typedef std::shared_ptr<const PeakDef> PeakConstPtr;
  typedef std::deque<PeakConstPtr> PeakDequeT;

  if( !background_spectrum || foreground_peaks.empty() )
    return background_auto_peaks;

  assert( fitPrefs );
  const PeakFitUtils::CoarseResolutionType prefs_det_type = fitPrefs
    ? fitPrefs->m_det_type
    : PeakFitUtils::CoarseResolutionType::Unknown;
  const PeakFitUtils::CoarseResolutionType det_type
    = (prefs_det_type != PeakFitUtils::CoarseResolutionType::Unknown)
        ? prefs_det_type
        : PeakFitUtils::coarse_det_type( background_spectrum, nullptr );

  const size_t nchannel = background_spectrum->num_gamma_channels();
  if( nchannel < 16 )
    return background_auto_peaks;

  const double min_energy = background_spectrum->gamma_channel_lower( 0 );
  const double max_energy = background_spectrum->gamma_channel_upper( nchannel - 1 );

  // Running background peak set we will augment (starts as a copy of the input).
  vector<PeakConstPtr> bg_peaks;
  if( background_auto_peaks )
    bg_peaks.insert( bg_peaks.end(), background_auto_peaks->begin(), background_auto_peaks->end() );

  // --- Step 1: collect unmatched, non-co-located foreground candidate energies, building a rough
  //     seed peak for each on the BACKGROUND (used both for causal clustering and as a fit seed). ---
  vector<PeakConstPtr> seeds;

  for( const PeakConstPtr &fg : foreground_peaks )
  {
    if( !fg || !fg->gausPeak() )
      continue;

    const double energy = fg->mean();
    const double fwhm = fg->fwhm();
    const double sigma = fg->sigma();
    if( (fwhm <= 0.0) || (sigma <= 0.0) || (energy <= min_energy) || (energy >= max_energy) )
      continue;

    // Skip if a background auto peak already matches (same rule the elevation test uses).
    bool matched = false;
    for( const PeakConstPtr &bg : bg_peaks )
    {
      const double ediff = fabs( energy - bg->mean() );
      const double tol = 0.75 * 0.5 * (fwhm + bg->fwhm());
      if( ediff < tol )
      {
        matched = true;
        break;
      }
    }//for( const PeakConstPtr &bg : bg_peaks )

    if( matched )
      continue;

    // Skip if we already have a seed essentially on top of this one.
    bool duplicate = false;
    for( const PeakConstPtr &s : seeds )
    {
      if( fabs( s->mean() - energy ) < 0.75*std::max( fwhm, s->fwhm() ) )
      {
        duplicate = true;
        break;
      }
    }//for( const PeakConstPtr &s : seeds )

    if( duplicate )
      continue;

    // Rough starting amplitude from the background gross counts over ~+-FWHM.
    const double gross = background_spectrum->gamma_integral( static_cast<float>( energy - fwhm ),
                                                              static_cast<float>( energy + fwhm ) );
    seeds.push_back( make_shared<PeakDef>( energy, sigma, std::max( 1.0, 0.5*gross ) ) );
  }//for( const PeakConstPtr &fg : foreground_peaks )

  if( seeds.empty() )
    return background_auto_peaks;

  // --- Step 2: cluster seeds using the SAME causal grouping the whole-spectrum search uses (7.5
  //     sigma), so a close doublet is reconstructed as ONE multi-peak ROI instead of two independent
  //     single-peak fits (one of which may drop the weaker line). ---
  const vector<vector<PeakConstPtr>> clusters = causilyDisconnectedPeaks( 7.5, true, seeds );

  // --- Step 3+4: for each cluster with a real background excess (cheap Currie pre-screen over the
  //     whole cluster span - side-bands OUTSIDE it, so a neighboring line doesnt bias the continuum
  //     and falsely reject a real one), fit ALL its seeds together on the background as a single
  //     multi-peak ROI, then merge into the running set (removing any superseded existing peaks so
  //     the ROI stays consistent). ---
  vector<PeakConstPtr> running = bg_peaks;   // running background peak set (starts as the base set)
  bool changed = false;

  for( const vector<PeakConstPtr> &cluster : clusters )
  {
    if( is_search_canceled( cancel_flag ) )
      break;
    if( cluster.empty() )
      continue;

    double span_lo = std::numeric_limits<double>::infinity();
    double span_hi = -std::numeric_limits<double>::infinity();
    for( const PeakConstPtr &s : cluster )
    {
      span_lo = std::min( span_lo, s->mean() - 1.25*s->fwhm() );
      span_hi = std::max( span_hi, s->mean() + 1.25*s->fwhm() );
    }

    if( (span_lo <= min_energy) || (span_hi >= max_energy) )
      continue;

    try
    {
      // Cheap Currie pre-screen over the whole cluster span.
      DetectionLimitCalc::CurrieMdaInput input;
      input.spectrum = background_spectrum;
      input.gamma_energy = static_cast<float>( 0.5*(span_lo + span_hi) );
      input.roi_lower_energy = static_cast<float>( span_lo );
      input.roi_upper_energy = static_cast<float>( span_hi );
      input.num_lower_side_channels = 4;
      input.num_upper_side_channels = 4;
      input.detection_probability = 0.95;
      input.additional_uncertainty = 0.0f;

      const DetectionLimitCalc::CurrieMdaResult curie = DetectionLimitCalc::currie_mda_calc( input );
      if( curie.source_counts <= curie.decision_threshold )
        continue;  // no real background excess in this region

      // Existing background peaks within the span are co-fit (and superseded) so the ROI stays whole.
      vector<PeakConstPtr> existing_in_span;
      for( const PeakConstPtr &p : running )
      {
        if( p && p->gausPeak() && (p->mean() >= span_lo) && (p->mean() <= span_hi) )
          existing_in_span.push_back( p );
      }

      // One shared continuum for the ROI, estimated from the background side channels.
      const shared_ptr<PeakContinuum> cont = make_shared<PeakContinuum>();
      cont->calc_linear_continuum_eqn( background_spectrum, 0.5*(span_lo + span_hi), span_lo, span_hi, 4, 4 );

      vector<PeakDef> fit_input;
      for( const PeakConstPtr &p : existing_in_span )
      {
        PeakDef s = *p;
        s.setContinuum( cont );
        s.setFitFor( PeakDef::Mean, true );
        s.setFitFor( PeakDef::Sigma, true );
        s.setFitFor( PeakDef::GaussAmplitude, true );
        fit_input.push_back( s );
      }

      for( const PeakConstPtr &seed : cluster )
      {
        // Skip a seed that essentially coincides with an existing peak already added above.
        bool coincides = false;
        for( const PeakConstPtr &p : existing_in_span )
        {
          if( fabs( p->mean() - seed->mean() ) < 0.5*seed->fwhm() )
          {
            coincides = true;
            break;
          }
        }
        if( coincides )
          continue;

        PeakDef s = *seed;
        s.setContinuum( cont );
        s.setFitFor( PeakDef::Mean, true );
        s.setFitFor( PeakDef::Sigma, true );
        s.setFitFor( PeakDef::GaussAmplitude, true );
        fit_input.push_back( s );
      }//for( const PeakConstPtr &seed : cluster )

      if( fit_input.empty() )
        continue;

      // Fit the whole ROI on the background as one group (large ncausality keeps the cluster
      //  together); peaks below the significance threshold are dropped.  Because the Currie
      //  pre-screen already confirmed a real excess, a modest significance threshold is used.
      const double ncausality = 7.5;
      const double stat_threshold = 2.0;       // minimum area significance (LM path)
      const double hypothesis_threshold = -1.0;
      const vector<PeakDef> fit = fitPeaksInRange( span_lo, span_hi, ncausality, stat_threshold,
                                        hypothesis_threshold, fit_input, background_spectrum,
                                        Wt::WFlags<PeakFitLM::PeakFitLMOptions>(), det_type );

      // Require netting at least one new peak beyond what was already there.
      if( fit.size() <= existing_in_span.size() )
        continue;

      // Remove the superseded existing peaks, then add the refit ROI peaks (which share `cont`).
      for( const PeakConstPtr &p : existing_in_span )
      {
        const vector<PeakConstPtr>::iterator pos = std::find( running.begin(), running.end(), p );
        if( pos != running.end() )
        {
          running.erase( pos );
          changed = true;
        }
      }

      for( const PeakDef &p : fit )
      {
        running.push_back( make_shared<PeakDef>( p ) );
        changed = true;
      }
    }catch( std::exception & )
    {
      // ROI near the spectrum edge or a fit error - just skip this cluster.
    }
  }//for( const vector<PeakConstPtr> &cluster : clusters )

  if( !changed )
    return background_auto_peaks;

  std::sort( running.begin(), running.end(), &PeakDef::lessThanByMeanShrdPtr );

  return std::make_shared<PeakDequeT>( running.begin(), running.end() );
}//recover_background_peaks_under_foreground(...)

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


void SavitzyGolayCoeffs::smooth_with_variance( const std::vector<float> &input,
                                               std::vector<float> &output,
                                               std::vector<float> &variance ) const
{
  const int nSamples = static_cast<int>( input.size() );
  const int nCoeffs = static_cast<int>( coeffs.size() );

  if( nSamples < nCoeffs || nSamples == 0 )
    throw std::runtime_error( "SavitzyGolayCoeffs::smooth_with_variance(...)\n\tInvalid input size" );

  output.clear();
  output.resize( nSamples );
  variance.clear();
  variance.resize( nSamples );

  for( int pos = 0; pos < nSamples; ++pos )
  {
    double sum = 0.0;
    double var = 0.0;

    for( int coeff = 0; coeff < nCoeffs; ++coeff )
    {
      int dataInd = pos - num_left + coeff;
      if( dataInd < 0 )
        dataInd = 0;
      else if( dataInd >= nSamples )
        dataInd = nSamples - 1;

      const double c = coeffs[coeff];
      const float N = input[dataInd];

      // Smoothed value
      sum += c * N;

      // Variance propagation for Poisson data: σ²(N) = N
      // For linear filter: σ²(output) = Σ c² × σ²(input) = Σ c² × N
      var += c * c * N;
    }

    output[pos] = sum;
    variance[pos] = var;
  }
}

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
  // Use group_peaks_by_roi for deterministic ordering (sorted by lowerEnergy)
  std::vector<std::pair<const PeakContinuum *, std::vector<PeakDef>>> groups
    = group_peaks_by_roi( input_peaks );

  for( auto &group : groups )
  {
    group.second[0].makeUniqueNewContinuum();
    std::shared_ptr<PeakContinuum> newcont = group.second[0].continuum();
    for( size_t i = 1; i < group.second.size(); ++i )
      group.second[i].setContinuum( newcont );
  }

  input_peaks.clear();
  for( const auto &group : groups )
  {
    for( const PeakDef &p : group.second )
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
                          std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                          vector<std::shared_ptr<PeakDef> > &answer,
                          double &chi2 )
{
  assert( fitPrefs );
  const PeakFitUtils::CoarseResolutionType det_type
    = PeakFitUtils::effective_det_type( fitPrefs, dataH, nullptr );
  const bool isHPGe = (det_type == PeakFitUtils::CoarseResolutionType::High);

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
  
  
  const size_t nSideBinsToAverage = 3;
  PeakContinuum::OffsetType offsetType = PeakContinuum::Linear;
  if( intputSharesContinuum && !answer.empty() )
    offsetType = answer[0]->continuum()->type();
  else if( nPeaks > 3 && (end_channel - start_channel) > 20 )  //20 is a WAG
    offsetType = PeakContinuum::Quadratic;
  
  double p0, p1;
  PeakContinuum::eqn_from_offsets( start_channel, end_channel, start_range,
                                  dataH, nSideBinsToAverage, nSideBinsToAverage, p1, p0 );
  
  
  double conteqn[6] = { p0, p1, 0.0, 0.0, 0.0, 0.0 };
  double initial_cont_area = PeakContinuum::offset_eqn_integral( conteqn, offsetType, start_range, end_range, start_range );
  
  const double totalpeakarea = (areaarea > initial_cont_area)
                                 ? (areaarea - initial_cont_area) : areaarea;
  
  
  vector<std::shared_ptr<const PeakDef>> inpeaks;
  std::shared_ptr<PeakContinuum> shared_continuum;
  switch( method )
  {
    case UniformInitialGuess:
    {
      for( int i = 0; i < nPeaks; ++i )
      {
        auto peak = make_shared<PeakDef>();
        
        const double amp = totalpeakarea / nPeaks;
        double mean = x0 + (i+0.5)*(x1-x0)/nPeaks;
        double sigma = 0.25*fabs(x1-x0) / nPeaks;
        
        if( detector && detector->hasResolutionInfo() )
          sigma = detector->peakResolutionSigma( mean );
        
        peak->setAmplitude( amp );
        peak->setMean( mean );
        peak->setSigma( sigma );
        
        if( i == 0 )
          shared_continuum = peak->continuum();
        else
          peak->setContinuum( shared_continuum );
        
        inpeaks.push_back( peak );
      }//for( int i = 0; i < nPeaks; ++i )
      
      break;
    }//case UniformInitialGuess:
      
    case FromDataInitialGuess:
    {
      intputSharesContinuum = false;
      std::vector< std::tuple<float,float,float> > candidates;  //{mean,sigma,area}
      secondDerivativePeakCanidates( dataH, fitPrefs, start_channel, end_channel, candidates );
      std::sort( begin(candidates), end(candidates),
                []( const tuple<float,float,float> &lhs, const tuple<float,float,float> &rhs) -> bool {
        return std::get<2>(lhs) > std::get<2>(rhs);
      } );
      
      double avrg_width = 0.0, avrg_amp = 0.0;
      for( auto i = begin(candidates); int(inpeaks.size()) < nPeaks && i != end(candidates); ++i )
      {
        if( get<0>(*i) >= x0 && get<0>(*i) <= x1 )
        {
          avrg_width += get<1>(*i);
          avrg_amp += get<2>(*i);
          
          auto peak = make_shared<PeakDef>(get<0>(*i), get<1>(*i), get<2>(*i));
          
          if( !shared_continuum )
            shared_continuum = peak->continuum();
          else
            peak->setContinuum( shared_continuum );
          
          inpeaks.push_back( peak );
        }//if( i->second.mean() >= x0 && i->second.mean() <= x1 )
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
        auto peak = make_shared<PeakDef>();
        const double amp = totalpeakarea / nPeaks;
        const double mean = x0 + (x1-x0)/(1+inpeaks.size());
        double sigma = avrg_width;
        if( detector && detector->hasResolutionInfo() )
          sigma = detector->peakResolutionSigma( mean );
        
        peak->setMean( mean );
        peak->setSigma( sigma );
        peak->setAmplitude( amp );
        
        if( i == 0 )
          shared_continuum = peak->continuum();
        else
          peak->setContinuum( shared_continuum );
        
        inpeaks.push_back( peak );
      }//for( int i = int(inpeaks.size()); i < nPeaks; ++i )
      
      break;
    }//case FromDataInitialGuess:
      
    case FromInputPeaks:
    {
      if( answer.size() != static_cast<size_t>(nPeaks) )
        throw runtime_error( "findPeaksInUserRange: invalid input for method=FromInputPeaks" );
      
      for( size_t i = 0; i < answer.size(); ++i )
      {
        auto peak = make_shared<PeakDef>( *answer[i] );
        if( i == 0 )
          shared_continuum = std::make_shared<PeakContinuum>( *peak->continuum() );
        peak->setContinuum( shared_continuum );
        inpeaks.push_back( peak );
      }
      
      answer.clear();
      break;
    }//case FromInputPeaks:
  }//switch( method )
  
  assert( shared_continuum );
  
  shared_continuum->setRange( x0, x1 );
  shared_continuum->setType( offsetType );
  
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> lm_fit_options;
  const double stat_threshold = 0.0;
  const double hypothesis_threshold = 0.0;
  
  vector<shared_ptr<const PeakDef>> results;
  
  PeakFitLM::fit_peaks_LM( results, inpeaks, dataH, stat_threshold, hypothesis_threshold, lm_fit_options, det_type );
  
  answer.clear();
  if( static_cast<int>(results.size()) == nPeaks )
  {
    for( const shared_ptr<const PeakDef> &p : results )
      answer.push_back( make_shared<PeakDef>(*p) );
  }else
  {
    cerr << "findPeaksInUserRange: fit " << results.size()
    << "peaks, but was requested " << nPeaks << " so will not return any peaks." << endl;
    
  }//if( two peaks are too close, or just really couldnt fit the peaks )
  
  if( !answer.empty() )
    chi2 = answer.front()->chi2dof();
  
  return;
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
                                     const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options,
                                     const PeakFitUtils::CoarseResolutionType det_type )
{
  //20120309: For the Ba133 example spectrum with default settings on my newer
  //          mac book pro, this function takes:
  //Single Thread:                   0.069735s wall, 0.070000s user + 0.000000s system = 0.070000s CPU (100.4%)
  //Multithreaded (phys cores only): 0.022624s wall, 0.080000s user + 0.000000s system = 0.080000s CPU (353.6%)
  //Multithreaded (logical cores)  : 0.019906s wall, 0.100000s user + 0.000000s system = 0.100000s CPU (502.3%)

  if( !data || (x1<x0) )
    return input_peaks;

  vector<shared_ptr<const PeakDef>> peak_from_LM;
  try
  {
    vector<shared_ptr<const PeakDef>> lm_input_peaks;
    for( const PeakDef &p : input_peaks )
      lm_input_peaks.push_back( make_shared<PeakDef>(p) );

    peak_from_LM = PeakFitLM::fit_peaks_in_range_LM( x0, x1, ncausality, stat_threshold, hypothesis_threshold,
                                                    lm_input_peaks, data, fit_options, det_type );
  }catch( std::exception &e )
  {
    // A Minuit2 based fit used to be the fallback here; with it gone, leave the peaks as they
    //  were, which is what this function also does for invalid input.
    cerr << "PeakFitLM::fit_peaks_in_range_LM threw: '" << e.what() << "'." << endl;
    assert( 0 );
    return input_peaks;
  }//try / catch

  input_peaks.clear();
  for( const shared_ptr<const PeakDef> &p : peak_from_LM )
    input_peaks.push_back( *p );

  //  cout << "Fit took: " << timer.format() << endl;

  return input_peaks;
}//std::vector<PeakDef> fitPeaksInRange(...)


//Defined below; SNIP clipping with a per-channel final clip half-width.
static const char *calculate_continuum_variable_window( float *spectrum, const int ssize,
                                                        const std::vector<int> &channel_windows,
                                                        const int filter_order,
                                                        const int presmooth_halfwidth,
                                                        const bool lls,
                                                        const int clip_begin,
                                                        const int clip_end );

std::shared_ptr<Measurement> estimateContinuum( std::shared_ptr<const Measurement> data )
{
  if( !data )
    throw runtime_error( "estimateContinuum: invalid data" );

  const size_t nchannels = data->num_gamma_channels();

  auto source = std::make_shared<vector<float>>( nchannels );
  for( size_t i = 0; i < nchannels; ++i )
    (*source)[i] = data->gamma_channel_content(i);

  // The fixed clip half-window this app has historically used (was `numIteration` of the
  // since-removed ROOT-derived calculateContinuum; the clean-room kernel below reproduces it to
  // float rounding).  Requires 2*window+1 <= nchannels; for very small spectra there is nothing
  // to clip, so - matching the old wrapper's silent no-op there - return the data unchanged.
  if( nchannels >= 5 )
  {
    const int window = std::min( 125, (static_cast<int>(nchannels) - 1) / 2 );
    const std::vector<int> windows( nchannels, window );
    const char *err = calculate_continuum_variable_window( &(source->operator[](0)),
                          static_cast<int>(nchannels), windows, 6, 0, false,
                          0, static_cast<int>(nchannels) );
    assert( !err );
    (void)err;
  }//if( nchannels >= 5 )

  auto background = make_shared<Measurement>();
  *background = *data;

  assert( background->num_gamma_channels() == nchannels );

  background->set_gamma_counts( source, data->live_time(), data->real_time() );

  return background;
}//std::shared_ptr<Measurement> estimateContinuum( std::shared_ptr<const Measurement> data )


/** SNIP continuum estimator (Statistics-sensitive Non-linear Iterative Peak-clipping),
 implemented from the published algorithm descriptions - C.G. Ryan et al., Nucl. Instr. Meth.
 B 34 (1988) 396-402, and M. Morhac et al., Nucl. Instr. Meth. A 401 (1997) 113-132 - with an
 increasing clip window and the higher-order (up to 6th) clipping filters, generalized here to
 a per-channel final clip half-width.

 The algorithm: iteratively replace each channel's value v(j) with
 min( v(j), P_w(j) ), where P_w(j) estimates what a smooth continuum would be at j from
 symmetric samples a clip half-width w away.  Peaks (narrower than w) are eaten away from the
 top over the iterations; the smooth continuum under them survives.  The continuum estimates
 come from requiring the (2p)-th central finite difference of a smooth curve to vanish, then
 solving for the center value:
   order 2:  E2 = ( v(j-w) + v(j+w) ) / 2
   order 4:  E4 = ( 4*(v(j-h2)+v(j+h2)) - (v(j-2*h2)+v(j+2*h2)) ) / 6,      h2 = w/2
   order 6:  E6 = ( 15*(v(j-h3)+v(j+h3)) - 6*(v(j-2*h3)+v(j+2*h3))
                    + (v(j-3*h3)+v(j+3*h3)) ) / 20,                          h3 = w/3
 (the coefficients are the binomial coefficients of the vanishing central difference).
 `filter_order` selects which estimates enter the composite clip, max(E2[,E4[,E6]]):
   - order 2 clips against E2 alone.  The taps sit at exactly +/-w, so with w sized to the peak
     half-width the peak is erased cleanly; the tradeoff is that E2 rounds sharp continuum
     corners (Compton edges) toward the two-sided average.
   - order 4/6 add E4/E6, whose taps sit at fractional distances (w/2, and w/3, 2w/3).  Those
     follow genuine continuum curvature, but for a window sized to ~1-2 FWHM the fractional taps
     land INSIDE a peak of that width and pull max(...) up, leaving a residual under strong
     peaks.  Higher orders therefore only pay off when w is much larger than the peak width
     (e.g. the legacy fixed 125-channel window).  Note E4/E6 also degenerate to v(j) for
     w < 2 / w < 3 (h2/h3 truncate to 0), so an order 4/6 filter does not clip until w >= 2/3.

 `channel_windows[j]` is the final clip half-width, in channels, for channel j.  Each channel's
 clip distance grows by +1 per iteration and saturates at its own window:
 d_j(m) = min( m, channel_windows[j] ) for iteration m = 1..max(channel_windows) (the standard
 increasing-window schedule, applied per channel).  A channel is only clipped by windows that
 fully fit inside the spectrum.

 If `presmooth` is true, a 3-channel boxcar is applied to the working spectrum before clipping,
 so the min-filter does not lock onto the downward noise fluctuations.

 If `lls` is true the clipping is performed in log-log-sqrt space (the LLS operator of Ryan et
 al., v = ln(ln(sqrt(y+1)+1)+1)), which compresses the dynamic range so high-count noise biases
 the min-filter less; the result is transformed back to counts.

 Returns nullptr on success, else an error message.
 */
static const char *calculate_continuum_variable_window( float *spectrum, const int ssize,
                                                        const std::vector<int> &channel_windows,
                                                        const int filter_order,
                                                        const int presmooth_halfwidth,
                                                        const bool lls,
                                                        const int clip_begin,
                                                        const int clip_end )
{
  if( ssize <= 0 )
    return "Wrong Parameters";

  if( (filter_order != 2) && (filter_order != 4) && (filter_order != 6) )
    return "filter_order must be 2, 4, or 6";

  if( static_cast<int>(channel_windows.size()) != ssize )
    return "channel_windows size must match spectrum size";

  // Only channels in [clip_begin, clip_end) are clipped; the rest are left equal to the input
  // data.  Taps never read outside this range, so a huge feature just outside it (e.g. the
  // detector turn-on below the lower spectroscopic extent) cannot pull the in-range continuum up.
  const int cbeg = std::max( 0, clip_begin );
  const int cend = std::min( ssize, clip_end );
  if( cbeg >= cend )
    return nullptr;  // nothing to clip; spectrum left unchanged

  const int span = cend - cbeg;

  int max_window = 0;
  for( int j = cbeg; j < cend; ++j )
  {
    const int w = channel_windows[j];
    if( w < 1 )
      return "Width of Clipping Window Must Be Positive";
    max_window = std::max( max_window, w );
  }

  if( span < (2*max_window + 1) )
    return "Too Large Clipping Window";

  // Double buffer: `prev` holds the previous iteration's values the clip filter reads;
  // `cur` receives this iteration's clipped values.
  std::vector<float> prev_buf( static_cast<size_t>(ssize) ), cur_buf( static_cast<size_t>(ssize) );
  float *prev = prev_buf.data(), *cur = cur_buf.data();

  for( int i = 0; i < ssize; i++ )
    prev[i] = cur[i] = spectrum[i];

  if( presmooth_halfwidth > 0 )
  {
    // (2*presmooth_halfwidth+1)-channel boxcar over the clip range, so the min-filter locks onto
    // the local mean rather than the downward Poisson fluctuations (helps low-count regions).
    for( int j = cbeg; j < cend; j++ )
    {
      float sum = 0.0f, num = 0.0f;
      for( int w = std::max(j - presmooth_halfwidth, cbeg);
           w <= std::min(j + presmooth_halfwidth, cend - 1); w++ )
      {
        sum += spectrum[w];
        num += 1.0f;
      }
      prev[j] = cur[j] = sum / num;
    }//for( int j = cbeg; j < cend; j++ )
  }//if( presmooth_halfwidth > 0 )

  if( lls )
  {
    // LLS operator of Ryan et al.: compresses dynamic range so high-count noise biases the
    // min-filter less; inverted after clipping.
    for( int j = cbeg; j < cend; j++ )
    {
      const float y = std::max( prev[j], 0.0f );
      prev[j] = cur[j] = std::log( std::log( std::sqrt(y + 1.0f) + 1.0f ) + 1.0f );
    }
  }//if( lls )

  for( int m = 1; m <= max_window; m++ )
  {
    for( int j = cbeg; j < cend; j++ )
    {
      const int w = std::min( m, channel_windows[j] );
      if( ((j - w) < cbeg) || ((j + w) >= cend) )
        continue;

      // Continuum estimates at j from symmetric samples, per the vanishing-central-difference
      // filters described in the function comment.
      float est = 0.5f * ( prev[j - w] + prev[j + w] );  // E2

      if( filter_order >= 4 )
      {
        const int h2 = w / 2;
        const float e4 = ( 4.0f * (prev[j - h2] + prev[j + h2])
                           - (prev[j - 2*h2] + prev[j + 2*h2]) ) / 6.0f;
        est = std::max( est, e4 );
      }

      if( filter_order >= 6 )
      {
        const int h3 = w / 3;
        const float e6 = ( 15.0f * (prev[j - h3] + prev[j + h3])
                           - 6.0f * (prev[j - 2*h3] + prev[j + 2*h3])
                           + (prev[j - 3*h3] + prev[j + 3*h3]) ) / 20.0f;
        est = std::max( est, e6 );
      }

      cur[j] = std::min( prev[j], est );
    }//for( int j = cbeg; j < cend; j++ )

    std::copy( cur + cbeg, cur + cend, prev + cbeg );
  }//for( int m = 1; m <= max_window; m++ )

  if( lls )
  {
    for( int j = cbeg; j < cend; j++ )
    {
      const float ee = std::exp( std::exp( cur[j] ) - 1.0f ) - 1.0f;
      cur[j] = std::max( ee*ee - 1.0f, 0.0f );
    }
  }//if( lls )

  std::copy( cur + cbeg, cur + cend, spectrum + cbeg );

  return nullptr;
}//calculate_continuum_variable_window(...)


std::shared_ptr<Measurement> estimateContinuum( std::shared_ptr<const Measurement> data,
                                                const std::function<double(double)> &fwhm_at_energy,
                                                const double num_fwhm_window,
                                                const int filter_order,
                                                const int presmooth_halfwidth,
                                                const bool lls,
                                                const double restrict_lower_energy,
                                                const double restrict_upper_energy )
{
  if( !data )
    throw runtime_error( "estimateContinuum: invalid data" );

  if( !fwhm_at_energy )
    throw runtime_error( "estimateContinuum: invalid fwhm function" );

  if( !(num_fwhm_window > 0.0) )
    throw runtime_error( "estimateContinuum: num_fwhm_window must be positive" );

  const size_t nchannels = data->num_gamma_channels();
  if( nchannels < 16 )
    throw runtime_error( "estimateContinuum: too few channels" );

  // Convert num_fwhm_window x FWHM(E) to a per-channel half-window in channels.
  // Floor of 2 keeps E2's taps distinct from the center; the upper clamp keeps the largest
  // window valid for the clip loop's edge handling.
  const int min_window = 2;
  const int max_window = std::max( min_window, static_cast<int>(nchannels/2) - 1 );

  std::vector<int> windows( nchannels, 0 );  //0 marks "no valid FWHM here" until filled below
  for( size_t i = 0; i < nchannels; ++i )
  {
    const double lower = data->gamma_channel_lower( i );
    const double upper = data->gamma_channel_upper( i );
    const double width = upper - lower;
    const double fwhm = fwhm_at_energy( 0.5*(lower + upper) );

    if( (width > 0.0) && std::isfinite(fwhm) && (fwhm > 0.0) )
    {
      const double w = std::round( num_fwhm_window * fwhm / width );
      windows[i] = std::clamp( static_cast<int>(w), min_window, max_window );
    }
  }//for( size_t i = 0; i < nchannels; ++i )

  // Fill channels without a valid FWHM from their nearest valid neighbor.
  int last_valid = 0;
  for( size_t i = 0; i < nchannels; ++i )
  {
    if( windows[i] > 0 )
      last_valid = windows[i];
    else if( last_valid > 0 )
      windows[i] = last_valid;
  }

  last_valid = 0;
  for( size_t i = nchannels; i > 0; --i )
  {
    if( windows[i-1] > 0 )
      last_valid = windows[i-1];
    else if( last_valid > 0 )
      windows[i-1] = last_valid;
  }

  if( windows[0] <= 0 )
    throw runtime_error( "estimateContinuum: no channel had a valid FWHM" );

  // Optional restriction: clip the continuum only within [restrict_lower_energy,
  // restrict_upper_energy].  Channels outside are left equal to the data (so a caller that gates
  // on data-minus-continuum sees zero excess there), which keeps a huge out-of-range feature -
  // e.g. the detector turn-on below the lower spectroscopic extent - from pulling the in-range
  // continuum up through the min-filter's taps.
  int clip_begin = 0, clip_end = static_cast<int>( nchannels );
  if( (restrict_upper_energy > restrict_lower_energy) && (restrict_lower_energy > 0.0) )
  {
    clip_begin = static_cast<int>( data->find_gamma_channel( static_cast<float>(restrict_lower_energy) ) );
    clip_end = static_cast<int>( data->find_gamma_channel( static_cast<float>(restrict_upper_energy) ) ) + 1;
    clip_begin = std::clamp( clip_begin, 0, static_cast<int>(nchannels) );
    clip_end = std::clamp( clip_end, 0, static_cast<int>(nchannels) );
  }

  auto source = std::make_shared<vector<float>>( nchannels );
  for( size_t i = 0; i < nchannels; ++i )
    (*source)[i] = data->gamma_channel_content( i );

  const char *err = calculate_continuum_variable_window( &(source->operator[](0)),
                              static_cast<int>(nchannels), windows, filter_order,
                              presmooth_halfwidth, lls, clip_begin, clip_end );
  if( err )
    throw runtime_error( "estimateContinuum: " + std::string(err) );

  auto background = make_shared<Measurement>();
  *background = *data;

  assert( background->num_gamma_channels() == nchannels );

  background->set_gamma_counts( source, data->live_time(), data->real_time() );

  return background;
}//estimateContinuum( data, fwhm_at_energy, num_fwhm_window, presmooth )



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
        const double ncontinuum = continuum->offset_integral( xbinlow, xbinup, data, peaks );

        const double npeak = gauss_counts[i];
        const double uncert = ndata > PEAK_FIT_MIN_CHANNEL_UNCERT ? sqrt(ndata) : 1.0;
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
          const double ncontinuum = continuum->offset_integral( xbinlow, xbinup, data, peaks );
          predicted[channel] += ncontinuum + peak->gauss_integral( xbinlow, xbinup );
        }//for( int bin = xlowbin; bin <= xhighbin; ++bin )
      }//for( const PeakShrdPtr &peak : peaks )
      
      for( map<size_t,double>::iterator i = predicted.begin(); i != predicted.end(); ++i )
      {
        const double peakarea = i->second;
        const double ndata = data->gamma_channel_content( i->first );
        const double uncert = ndata > PEAK_FIT_MIN_CHANNEL_UNCERT ? sqrt(ndata) : 1.0;
        const double chi = (ndata-peakarea)/uncert;
        chi2 += chi*chi;
      }//for( int bin = minbin; bin <= maxbin; ++bin )
    }//if( continuum->energyRangeDefined() ) / else
  }//for( const ContToPeakMap_t::value_type &vt : contToPeakMap )
  
  return chi2;
}//chi2_for_region(...)



vector<shared_ptr<const PeakDef>> refitPeaksThatShareROI( const std::shared_ptr<const Measurement> &dataH,
                                   const std::shared_ptr<const DetectorPeakResponse> &detector,
                                   const PeakShrdVec &inpeaks,
                                   const PeakFitUtils::CoarseResolutionType det_type,
                                   const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options )
{
  return PeakFitLM::refitPeaksThatShareROI_LM( dataH, detector, inpeaks, det_type, fit_options );
}//PeakShrdVec refitPeaksThatShareROI(...)


double evaluate_chi2dof_for_range( const std::vector<PeakDef> &peaks,
                                  const std::shared_ptr<const Measurement> &dataH,
                                  const double startx,
                                  const double endx )
{
  const size_t lowerchannel = dataH->find_gamma_channel( startx );
  const size_t upperchannel = dataH->find_gamma_channel( endx );
  
  // Build a map of continuum to peaks for CDF step types
  map<const PeakContinuum *, vector<const PeakDef *>> cont_to_peaks;
  for( size_t j = 0; j < peaks.size(); ++j )
    cont_to_peaks[peaks[j].continuum().get()].push_back( &peaks[j] );

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
        const vector<const PeakDef *> &roi_peaks = cont_to_peaks[contptr];
        y_pred += contptr->offset_integral( x0, x1, dataH, roi_peaks.data(), roi_peaks.size() );
      }
    }//for( size_t j = 0; j < peaks.size(); ++j )
    
    if( y_pred < 0.0 )
      y_pred = 0.0;
    
    const double uncert = (y > PEAK_FIT_MIN_CHANNEL_UNCERT ? sqrt(y) : 1.0);
    chi2 += std::pow( (y_pred - y) / uncert, 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2 / (1 + upperchannel - lowerchannel);
}//double evaluate_chi2dof_for_range(...);


//refit_for_new_roi(): resultPeaks empty upon error.
void refit_for_new_roi( std::vector< std::shared_ptr<const PeakDef> > originalPeaks,
                       const std::shared_ptr<const Measurement> &dataH,
                       const double new_lower_roi,
                       const double new_roi_upper,
                       const PeakFitUtils::CoarseResolutionType det_type,
                       std::vector<PeakDef> &resultPeaks )
{
  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> lm_fit_options
                                  = PeakFitLM::PeakFitLMOptions::MediumRefinementOnly;
  const double stat_threshold = 0.0;
  const double hypothesis_threshold = 0.0;

  vector<shared_ptr<const PeakDef>> results;
  PeakFitLM::fit_peaks_LM( results, originalPeaks, dataH,
                          stat_threshold, hypothesis_threshold, lm_fit_options, det_type );
  
  resultPeaks.clear();
  
  if( !results.empty() && (results.size() == originalPeaks.size()) )
  {
    for( const shared_ptr<const PeakDef> &p : results )
      resultPeaks.push_back( *p );
  }
  
  return;
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
          PeakFitUtils::find_spectroscopic_extent( data, speclower, specupper );
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
                                 const PeakFitUtils::CoarseResolutionType det_type,
                                 const std::shared_ptr<const SpecUtils::Measurement> &meas,
                                 float &min_sigma_width_kev,
                                 float &max_sigma_width_kev )
{
  //For this function I took approximately the best and worst resolution
  //  detector response functions for lowres detectors in GADRAS, and this
  //  function just gives a mutliple of these values.
  //The arrays in this function were generated using the
  //  print_detector_sigma_range() function in developcode.cpp.
  const bool highres = (det_type == PeakFitUtils::CoarseResolutionType::High);
  const bool unknown = (det_type == PeakFitUtils::CoarseResolutionType::Unknown);

  const float max_width_multple = highres ? 4.0f : 3.0f;
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
  
  // Helper to interpolate from the lookup tables
  auto interp = [index, frac]( const float *arr ) -> float {
    return arr[index] + frac * (arr[index+1] - arr[index]);
  };

  const float sigma_to_fwhm = 2.3548f;

  switch( det_type )
  {
    case PeakFitUtils::CoarseResolutionType::High:
    {
      const float lowerres = interp( highres_smallest_expected_sigma );
      const float largerres = interp( highres_largest_expected_sigma );
      min_sigma_width_kev = min_width_multiple * lowerres;
      max_sigma_width_kev = max_width_multple * largerres;
      break;
    }

    case PeakFitUtils::CoarseResolutionType::Low:
    {
      const float lowerres = interp( lowres_smallest_expected_sigma );
      const float largerres = interp( lowres_largest_expected_sigma );
      min_sigma_width_kev = min_width_multiple * lowerres;
      max_sigma_width_kev = max_width_multple * largerres;
      break;
    }

    case PeakFitUtils::CoarseResolutionType::LaBr:
    {
      const float typical_sigma = PeakFitUtils::labr_fwhm_fcn( energy ) / sigma_to_fwhm;
      min_sigma_width_kev = 0.5f * typical_sigma;
      max_sigma_width_kev = 2.0f * typical_sigma;
      break;
    }

    case PeakFitUtils::CoarseResolutionType::CZT:
    {
      // "Good" CZT (M400): {2.7, 0.699, 0.753}
      static const std::vector<float> czt_good_coefs{ 2.7f, 0.699f, 0.753f };
      const float good_sigma = DetectorPeakResponse::peakResolutionFWHM( energy,
                    DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn,
                    czt_good_coefs ) / sigma_to_fwhm;
      // "General" CZT (Kromek GR1): {8.95, 2.39, 0.344}
      const float general_sigma = PeakFitUtils::czt_fwhm_fcn( energy ) / sigma_to_fwhm;

      min_sigma_width_kev = 0.5f * good_sigma;
      max_sigma_width_kev = 2.0f * general_sigma;
      break;
    }

    case PeakFitUtils::CoarseResolutionType::MedRes:
    {
      // MedRes could be either LaBr or CZT; use range spanning both
      const float labr_sigma = PeakFitUtils::labr_fwhm_fcn( energy ) / sigma_to_fwhm;

      static const std::vector<float> czt_good_coefs{ 2.7f, 0.699f, 0.753f };
      const float czt_good_sigma = DetectorPeakResponse::peakResolutionFWHM( energy,
                    DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn,
                    czt_good_coefs ) / sigma_to_fwhm;
      const float czt_general_sigma = PeakFitUtils::czt_fwhm_fcn( energy ) / sigma_to_fwhm;

      min_sigma_width_kev = 0.5f * std::min( labr_sigma, czt_good_sigma );
      max_sigma_width_kev = 2.0f * std::max( labr_sigma, czt_general_sigma );
      break;
    }

    case PeakFitUtils::CoarseResolutionType::LowOrMedRes:
    {
      // Could be NaI, CsI, LaBr, or CZT; use range from best CZT to worst NaI
      static const std::vector<float> czt_good_coefs{ 2.7f, 0.699f, 0.753f };
      const float czt_good_sigma = DetectorPeakResponse::peakResolutionFWHM( energy,
                    DetectorPeakResponse::ResolutionFnctForm::kGadrasResolutionFcn,
                    czt_good_coefs ) / sigma_to_fwhm;

      const float lowres_max = interp( lowres_largest_expected_sigma );

      min_sigma_width_kev = 0.5f * czt_good_sigma;
      max_sigma_width_kev = 3.0f * lowres_max;
      break;
    }

    case PeakFitUtils::CoarseResolutionType::Unknown:
    {
      // Use the widest possible range: min from HPGe, max from low-res
      const float hp_min_mult = (energy < 50.0f ? 0.25f : 0.35f);
      const float hp_lowerres = interp( highres_smallest_expected_sigma );
      const float lr_largerres = interp( lowres_largest_expected_sigma );

      min_sigma_width_kev = hp_min_mult * hp_lowerres;
      max_sigma_width_kev = 3.0f * lr_largerres;
      break;
    }
  }//switch( det_type )
  
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
                                 std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                                 const PeakShrdVec &inpeaks )
{
  assert( fitPrefs );
  const PeakFitUtils::CoarseResolutionType det_type
    = PeakFitUtils::effective_det_type( fitPrefs, dataH, nullptr );
  const bool isHPGe = (det_type == PeakFitUtils::CoarseResolutionType::High);

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
  expected_peak_width_limits( x, det_type, dataH, min_sigma_width_kev, max_sigma_width_kev );

  
  
  const size_t lower_reasonable_channel = dataH->find_gamma_channel( x - 20*max_sigma_width_kev );
  const size_t upper_reasonable_channel = dataH->find_gamma_channel( x + 20*max_sigma_width_kev );
  
  //cout << "For x=" << x << " keV, {lowchannel=" << lowchannel
  //     << ", lower_reasonable_channel=" << lower_reasonable_channel
  //     << "}, {" << highchannel << ", upper_reasonable_channel=" << upper_reasonable_channel << "}"
  //     << endl;
  
  lowchannel = std::max(lowchannel, lower_reasonable_channel);
  highchannel = std::min( highchannel, upper_reasonable_channel );
  
  
  const vector<PeakPtr> candidates
       = secondDerivativePeakCanidatesWithROI( dataH, fitPrefs, lowchannel, highchannel );
  

  float min_sigma, max_sigma;
  expected_peak_width_limits( x, det_type, dataH, min_sigma, max_sigma );

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
        
        vector<PeakDef> newfitpeaks;
        refit_for_new_roi( inpeaks, dataH, finalLowerE, startx, PeakFitUtils::CoarseResolutionType::LowOrMedRes, newfitpeaks );
        
        const double newChi2Dof = newfitpeaks.empty() ? DBL_MAX : evaluate_chi2dof_for_range( newfitpeaks, dataH, finalLowerE, startx );
        const double totalOldlen  = finalUpperE - finalLowerE;
        const double totalNewLen  = startx - finalLowerE;
        const double totalRemovedFrac = (totalOldlen - totalNewLen) / totalOldlen;
        double oldChi2Body = (chi2Dof - totalRemovedFrac*tailChi2Dof) /( 1.0 - totalRemovedFrac );
        
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cout) << "newChi2Dof=" << newChi2Dof << ", oldchi2dof=" << chi2Dof
             << ", oldChi2Body=" << oldChi2Body << "\n";
#endif
        
        //1.25 is arbitrary
        if( !newfitpeaks.empty() && (newChi2Dof < 1.25*oldChi2Body) )
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
        vector<PeakDef> newfitpeaks;
        refit_for_new_roi( inpeaks, dataH, endx, finalUpperE, PeakFitUtils::CoarseResolutionType::LowOrMedRes, newfitpeaks );
       
        const double newChi2Dof = newfitpeaks.empty() ? DBL_MAX : evaluate_chi2dof_for_range( newfitpeaks, dataH, endx, finalUpperE );
        const double totalOldlen  = finalUpperE - finalLowerE;
        const double totalNewLen  = finalUpperE - endx;
        const double totalRemovedFrac = (totalOldlen - totalNewLen) / totalOldlen;
        double oldChi2Body = (chi2Dof - totalRemovedFrac*tailChi2Dof) /( 1.0 - totalRemovedFrac );
          
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cout) << "newChi2Dof=" << newChi2Dof << ", oldchi2dof=" << chi2Dof << ", oldChi2Body=" << oldChi2Body << "\n";
#endif
        
        //1.25 is arbitrary
        if( !newfitpeaks.empty() && (newChi2Dof < 1.25*oldChi2Body) )
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
        
        vector<PeakDef> newfitpeaks;
          
        PeakShrdVec peaksToReFit;
        for( const PeakDef &p : fitpeaks )
          peaksToReFit.push_back( std::make_shared<PeakDef>(p) );
          
        const vector<double> originalFitPars, originalFitErrors;
          
        refit_for_new_roi( peaksToReFit, dataH, finalLowerE, finalUpperE, PeakFitUtils::CoarseResolutionType::LowOrMedRes, newfitpeaks );
          
        const double newChi2Dof = newfitpeaks.empty() ? DBL_MAX
                                           : evaluate_chi2dof_for_range( newfitpeaks, dataH, finalLowerE, finalUpperE );
        
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
      
      
      vector<PeakDef> newfitpeaks;
      refit_for_new_roi( inpeaks, dataH, test_lower_roi, test_upper_roi, PeakFitUtils::CoarseResolutionType::High, newfitpeaks );
      const double origChi2ForNewRange = evaluate_chi2dof_for_range( fitpeaks, dataH, test_lower_roi, test_upper_roi );
      const double testChi2Dof = newfitpeaks.empty() ? DBL_MAX
                                     : evaluate_chi2dof_for_range( newfitpeaks, dataH, test_lower_roi, test_upper_roi );
      
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
        vector<PeakDef> newfitpeaks;
        
        PeakShrdVec peaksToReFit;
        for( const PeakDef &p : fitpeaks )
          peaksToReFit.push_back( std::make_shared<PeakDef>(p) );
        
        const vector<double> originalFitPars, originalFitErrors;
        
        refit_for_new_roi( peaksToReFit, dataH, finalLowerE, finalUpperE, PeakFitUtils::CoarseResolutionType::High, newfitpeaks );
        
        const double newChi2Dof = newfitpeaks.empty() ? DBL_MAX
                                            : evaluate_chi2dof_for_range( newfitpeaks, dataH, finalLowerE, finalUpperE );
        
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
    expected_peak_width_limits( mean, PeakFitUtils::CoarseResolutionType::LowOrMedRes, dataH, min_sigma, max_sigma );

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
      const double uncert = (y > PEAK_FIT_MIN_CHANNEL_UNCERT ? std::sqrt(y) : 1.0);
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
    // This function is only called for single-peak ROIs, so this peak is the only ROI peer.
    const PeakDef *peak_ptr = peak.get();
    const double conttipval = peak->continuum()->offset_integral( mean-0.1*sigma, mean+0.1*sigma, dataH, &peak_ptr, 1 );
    const double cont2val = peak->continuum()->offset_integral( mean+1.9*sigma, mean+2.1*sigma, dataH, &peak_ptr, 1 );
    const double contdiff = conttipval - cont2val;
    const double peakdiff = ptipval - p2val;

    const double max_relative_continuum_slope = 1.25;
    const double minuncert_apply_cont_slope_test = 0.05;

    const double roi_lower = peak->lowerX();
    const double roi_upper = peak->upperX();

    const double below_roi_cont_area = peak->continuum()->offset_integral( roi_lower - sigma, roi_lower, dataH, &peak_ptr, 1 );
    const double above_roi_cont_area = peak->continuum()->offset_integral( roi_upper, roi_upper + sigma, dataH, &peak_ptr, 1 );
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
      const double val = continuum->offset_integral( energies[i], energies[i+1], dataH, fitpeaks );
      if( val < minval )
      {
        minval = val;
        minchanel = i;
      }
    }//for( size_t i = startchannel; i <= endchannel; ++i )

    const double lowedgeval = continuum->offset_integral( energies[startchannel], energies[startchannel+1], dataH, fitpeaks );
    const double highedgeval = continuum->offset_integral( energies[endchannel], energies[endchannel+1], dataH, fitpeaks );
  
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
      const double uncert = (y > PEAK_FIT_MIN_CHANNEL_UNCERT ? sqrt(y) : 1.0);
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
    vector<PeakDef> withoutResultPeaks;
    const double withNewChi2Dof = newpeak->chi2dof();
    
    refit_for_new_roi( otherpeak, dataH, lx, ux, PeakFitUtils::CoarseResolutionType::LowOrMedRes, withoutResultPeaks );
    
    const double withoutNewChi2Dof = withoutResultPeaks.empty() ? DBL_MAX
                                                  : evaluate_chi2dof_for_range( withoutResultPeaks, dataH, lx, ux );
    
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
    expected_peak_width_limits( mean, PeakFitUtils::CoarseResolutionType::High, dataH, min_sigma, max_sigma );

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
  expected_peak_width_limits( (float)mean, PeakFitUtils::CoarseResolutionType::High, dataH, min_sigma, max_sigma );
  
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
    const bool hasextent = PeakFitUtils::find_spectroscopic_extent(
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
    
    
    
    const PeakDef *peak_ptr_c = peak.get();
    const double lcont = peak->continuum()->offset_integral( lxl, lxu, dataH, &peak_ptr_c, 1 );
    const double ucont = peak->continuum()->offset_integral( uxl, uxu, dataH, &peak_ptr_c, 1 );
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
      const double uncert = (y > PEAK_FIT_MIN_CHANNEL_UNCERT ? std::sqrt(y) : 1.0);
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
                                                       const std::shared_ptr<const std::deque<shared_ptr<const PeakDef>>> &auto_search_peaks,
                                                       std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                                                       std::shared_ptr<const std::atomic<bool>> cancel_flag )
{
  typedef std::shared_ptr<const PeakDef> PeakDefShrdPtr;

  assert( fitPrefs );
  const bool isHPGe = (PeakFitUtils::effective_det_type( fitPrefs, dataH, nullptr )
                       == PeakFitUtils::CoarseResolutionType::High);

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
                                              pixelPerKev, dataH, fitPrefs, inpeaks );
  
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

  // Protect against a degenerate ROI from the heuristics above (e.g., a candidate near the
  //  spectrum edge, of a low-statistics spectrum, can yield inverted limits, or limits that dont
  //  contain the candidate) - treat as no peak found, rather than feeding an invalid range to the
  //  peak fit (which will assert on it, in debug builds).
  if( std::isnan(roiLower) || std::isnan(roiUpper) || (roiUpper <= roiLower)
      || (mean0 < roiLower) || (mean0 > roiUpper) )
  {
#if( PERFORM_DEVELOPER_CHECKS )
    char buffer[256];
    snprintf( buffer, sizeof(buffer), "Got degenerate ROI [%f, %f] for candidate peak mean=%f,"
              " sigma=%f - skipping peak fit.", roiLower, roiUpper, mean0, sigma0 );
    log_developer_error( __func__, buffer );
#endif
    return pair<PeakShrdVec,PeakShrdVec>();
  }//if( ROI limits are degenerate )

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
  
  
  vector<double> lowerEnergies, upperEnergies;
  lowerEnergies.push_back( roiLower );
  upperEnergies.push_back( roiUpper );

  
  PeakShrdVec initialfitpeaks;
  PeakFitLM::fit_peak_for_user_click_LM( initialfitpeaks, dataH, coFitPeaks,
                             mean0, sigma0, area0, lowerEnergies[0], upperEnergies[0],
                             fitPrefs, drf, cancel_flag );
  
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
      

      
      initialfitpeaks.clear();
      PeakFitLM::fit_peak_for_user_click_LM( initialfitpeaks, dataH, coFitPeaks,
                                 mean0, sigma0, area0, lowerEnergies[0], upperEnergies[0],
                                 fitPrefs, drf, cancel_flag );
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
                                   std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                                   size_t start_channel,
                                   size_t end_channel,
                                   std::vector< std::tuple<float,float,float> > &results )
{
  assert( fitPrefs );
  const bool isHPGe = (PeakFitUtils::effective_det_type( fitPrefs, data, nullptr )
                       == PeakFitUtils::CoarseResolutionType::High);

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
                                                                            std::shared_ptr<const PeakFitDetPrefs> fitPrefs,
                                                          size_t start_channel,
                                                          size_t end_channel )
{
  assert( fitPrefs );
  const bool isHPGe = (PeakFitUtils::effective_det_type( fitPrefs, dataH, nullptr )
                       == PeakFitUtils::CoarseResolutionType::High);

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
      const PeakDef *peak_ptr_dbg = peak.get();
      double cont_area = peak->continuum()->offset_integral( lowerEnengy, upperEnergy, dataH, &peak_ptr_dbg, 1 );
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
                              const vector<PeakDef *> &peaks,
                              const vector<float> *channel_count_uncerts )
{
  dof = chi2 = 0.0;
  if( peaks.empty() || !data || !data->channel_energies() || data->channel_energies()->empty() )
    return;

  if( channel_count_uncerts )
  {
    const size_t num_channels = data->num_gamma_channels();
    if( channel_count_uncerts->size() != num_channels )
      throw std::runtime_error( "get_chi2_and_dof_for_roi: channel_count_uncerts size doesn't match data" );

#if( PERFORM_DEVELOPER_CHECKS )
    for( size_t i = 0; i < num_channels; ++i )
    {
      const float uncert = (*channel_count_uncerts)[i];
      if( uncert <= 0.0f || std::isnan(uncert) || std::isinf(uncert) )
      {
        char buffer[256];
        snprintf( buffer, sizeof(buffer),
                  "get_chi2_and_dof_for_roi: invalid uncertainty at channel %zu: %g",
                  i, uncert );
        throw std::runtime_error( buffer );
      }
    }
#endif
  }
  
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

  vector<const PeakDef *> roi_peak_ptrs( peaks.size() );
  for( size_t pi = 0; pi < peaks.size(); ++pi )
    roi_peak_ptrs[pi] = peaks[pi];

  for( size_t i = 0; i < numchannel; ++i )
  {
    const size_t channel = startchannel + i;
    const double xbinlow = energies[channel];
    const double xbinup  = energies[channel+1];
    double nfitpeak = gauss_counts[i];

    const double ndata = data->gamma_channel_content( channel );
    const double ncontinuim = continuum->offset_integral( xbinlow, xbinup, data, roi_peak_ptrs.data(), roi_peak_ptrs.size() );

    // datauncert is variance (σ²), not standard deviation
    const double datauncert = channel_count_uncerts
                              ? static_cast<double>((*channel_count_uncerts)[channel] * (*channel_count_uncerts)[channel])
                              : std::max( ndata, 1.0 );

#if( PERFORM_DEVELOPER_CHECKS )
    if( channel_count_uncerts )
    {
      assert( datauncert > 0.0 );
      assert( !std::isnan(datauncert) && !std::isinf(datauncert) );
    }
#endif

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

  // FWHMs across an ROI are tied via a linear relation sigma(E) = a + b*E by default, so the FWHM
  //  cost is at most 2 parameters per ROI (matches the LM default in
  //  PeakFitDiffCostFunction::roi_sigma_parameter_count()).  That function returns every fitted
  //  sigma under the AllPeakFwhmIndependent option; this site cannot see the option flags, so it
  //  always assumes the linked convention.
  const int num_sigmas_fit = std::min( nfitsigma, 2 );

  // Continuum parameters consume DOF whether they are fit by the optimizer directly or solved
  //  analytically via linear least squares.  `fitForParameter()` is sized `num_parameters(...)`,
  //  which already includes any peak-CDF step coefficients, so counting it is the whole story -
  //  adding the step coefficients again would count them twice.
  int num_fit_continuum_pars = 0;
  for( const bool fit : continuum->fitForParameter() )
    num_fit_continuum_pars += (fit ? 1 : 0);

  // The LM fitter's counterpart is PeakFitDiffCostFunction::dof_for_roi() (src/PeakFitLM.cpp), and
  //  the two do NOT presently agree.  KNOWN DISCREPANCY, pre-existing and deliberately left alone
  //  (2026-09): on its least-squares path `dof_for_roi()` subtracts only
  //  `num_cdf_step_pars(type)` - zero for every non-CDF type - so it never charges for the
  //  polynomial coefficients `PeakFit::fit_amp_and_offset_imp(...)` solves, and overstates DOF by
  //  `num_linear_fit_pars(type)`: 1 for Constant, 2 for Linear, 4 for Cubic.  The convention used
  //  *here* is the correct one - a parameter that is fit costs a degree of freedom whether Ceres
  //  or the linear solve fits it - so the fix belongs on the LM side.  It is not applied because
  //  `PeakDef::chi2dof()` also gates automated peak acceptance (`max_chi2dof_roi`,
  //  `lowres_max_chi2dof`, and the one-vs-two-peak rule `twoPeaksChi2 <= chi2dof() + 1.2` in this
  //  file), so correcting it changes which peaks a search keeps, not just what is reported.
  //
  // Measured 2026-09-12 over 1362 auto-fit peaks (26 spectra, mixed HPGe/NaI): correcting the LM
  //  side raises reported chi2/DOF by 5-7% (median ratio 1.069 and 1.048 on the two corpora, range
  //  1.008-1.167) and, for every one of those peaks, changes *only* the divisor - chi2 itself did
  //  not move, and no peak was gained or lost.  The rule above that could still bite is the
  //  one-vs-two-peak margin, which is absolute rather than relative.
  //
  // Three further divergences are deliberate rather than defects, and should not be "fixed" in
  //  isolation:
  //   - skew: shared skew parameters are amortized across every ROI of an LM fit, so there is no
  //     correct per-ROI count to make here.
  //   - this function floors `dof` at 1.0; the LM side instead refuses the fit in its constructor.
  //   - the ROI channel window below uses a +-1e-7 keV epsilon that the LM side does not, which
  //     matters only when a ROI edge lands exactly on a channel boundary.
  const double num_channels = 1.0 + endchannel - startchannel;
  dof = num_channels - nfitamp - nfitmean - num_sigmas_fit - num_fit_continuum_pars;

  if( dof < 1.0 )
    dof = 1.0;
}//get_chi2_and_dof_for_roi(...)



double set_chi2_dof( std::shared_ptr<const Measurement> data,
                    std::vector<PeakDef> &fitpeaks,
                    const size_t startpeak, const size_t npeaks,
                    const vector<float> *channel_count_uncerts )
{
  if( channel_count_uncerts && data )
  {
    const size_t num_channels = data->num_gamma_channels();
    if( channel_count_uncerts->size() != num_channels )
      throw std::runtime_error( "set_chi2_dof: channel_count_uncerts size doesn't match data" );

    // Uncertainty validation will be done in get_chi2_and_dof_for_roi
  }

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
    get_chi2_and_dof_for_roi( chi2, dof, data, peakptrs, channel_count_uncerts );

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


std::vector<std::shared_ptr<const PeakDef>> peaksInRange( const double lowx,
                                   const double highx,
                                   const double nsigma,
                                   const std::vector<std::shared_ptr<const PeakDef>> &inputs )
{
  vector<shared_ptr<const PeakDef>> answer;

  for( const shared_ptr<const PeakDef> &peak : inputs )
  {
    const double peakmin  = (peak->gausPeak() ? (peak->mean() - nsigma*peak->sigma()) : peak->lowerX());
    const double peakmax  = (peak->gausPeak() ? (peak->mean() + nsigma*peak->sigma()) : peak->upperX());

    if( (peakmin<=highx) && (peakmax>=lowx) )  //covers case 1 and 3
      answer.push_back( peak );
    else if( (peakmin<=lowx) && (peakmax>=lowx) ) //covers case 2
      answer.push_back( peak );
  }//for( const PeakDef &peak : inputs )

  std::sort( begin(answer), end(answer), &PeakDef::lessThanByMeanShrdPtr );

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
  //Implementation using Eigen SVD for numerical stability
  const int poly_terms = polynomial_order + 1;
  Eigen::MatrixX<double> A( nbin, poly_terms );
  Eigen::VectorX<double> b( nbin );
  
  for( size_t row = 0; row < nbin; ++row )
  {
    const double uncert = (data[row] > PEAK_FIT_MIN_CHANNEL_UNCERT ? sqrt( data[row] ) : 1.0);
    b(row) = (data[row] > 0.0 ? sqrt( data[row] ) : 0.0);
    for( int col = 0; col < poly_terms; ++col )
      A(row,col) = std::pow( double(x[row]), double(col)) / uncert;
  }//for( int col = 0; col < poly_terms; ++col )
  
#if( EIGEN_VERSION_AT_LEAST( 3, 4, 1 ) )
  const Eigen::JacobiSVD<Eigen::MatrixX<double>,Eigen::ComputeThinU | Eigen::ComputeThinV> svd(A);
#else
  const Eigen::BDCSVD<Eigen::MatrixX<double>> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV );
#endif
  
  const Eigen::VectorXd a = svd.solve(b);
  
  const Eigen::MatrixX<double> A_transpose = A.transpose();
  const Eigen::MatrixX<double> alpha = A_transpose * A;
  const Eigen::MatrixX<double> C = alpha.inverse();
  
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
    const double uncert = (data[bin] > PEAK_FIT_MIN_CHANNEL_UNCERT ? sqrt( data[bin] ) : 1.0);
    chi2 += std::pow( (y_pred - data[bin]) / uncert, 2.0 );
  }//for( int bin = 0; bin < nbin; ++bin )
  
  return chi2;
}//double fit_to_polynomial(...)
        
        
double fit_amp_and_offset( const float *x, const float *data, const size_t nbin,
                                  const PeakContinuum::OffsetType cont_type,
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
  double * const dummy_channel_counts = nullptr;

  return PeakFit::fit_amp_and_offset_imp( x, data, nullptr, nbin, cont_type,
                  nullptr, ref_energy, means, sigmas, fixedAmpPeaks, skew_type, skew_parameters,
                                         amplitudes, continuum_coeffs,
                                         amplitudes_uncerts, continuum_coeffs_uncerts,
                                         dummy_channel_counts );
}//double fit_amp_and_offset(...)
        
        
        
bool chi2_significance_test( const PeakDef &peak,
                             const double withoutPeakDSigma,
                             const double chi2ratioRequired,
                             std::vector<PeakDef> other_peaks,
                             const std::shared_ptr<const Measurement> &data )
{
  if( !peak.gausPeak() )
    throw logic_error( "chi2_significance_test: can only evaluate gaussian defined peaks" );

  for( const PeakDef &p : other_peaks )
  {
    if( !p.gausPeak() )
      throw logic_error( "chi2_significance_test: all peaks in ROI must be gaussian defined" );
  }

  if( !data || (data->num_gamma_channels() < 7) || !data->energy_calibration() || !data->energy_calibration()->valid() )
    throw runtime_error( "chi2_significance_test: invalid spectrum" );


  std::shared_ptr<const PeakContinuum> cont = peak.continuum();
  const PeakContinuum::OffsetType offset_type = cont->type();

  // We will limit the test range to +-2.5 FWHM of the peak mean to do the Chi2 test so we capture,
  //  only the relevant part of the peak.  But the step-continuum peaks we need to keep the original
  //  energy range, or else it will alter each channels continuum values (and we will fail peak fits
  //  we shouldnt).
  //
  //TODO 20210726: I dont really remember why I limited this check to +-2.5 FWHM from peak mean, but
  //  I'm guessing its to avoid failing long continuums where the tails are badly matched to the
  //  data, but user wants in order to force things.  Should investigate/document more thoroughly.
  const bool is_step = PeakContinuum::is_step_continuum( offset_type );

  //For section of code below here, please see notes in searchForPeaks()
  other_peaks.erase( remove_if( begin(other_peaks), end(other_peaks), [&peak](const PeakDef &other) {
          return fabs(other.mean() - peak.mean()) < 0.0001 && fabs(other.sigma() - peak.sigma()) < 0.0001;
   } ), end(other_peaks) );

  //Some basic quality checks
  if( IsNan(peak.mean()) || IsInf(peak.mean()) )
    return false;

  if( IsNan(peak.sigma()) || IsInf(peak.sigma()) )
    return false;

  if( IsNan(peak.amplitude()) || IsInf(peak.amplitude()) )
    return false;

  if( peak.fitFor(PeakDef::CoefficientType::GaussAmplitude) && ((peak.amplitude() <= 0.0)) )
    return false;

  if( peak.fitFor(PeakDef::CoefficientType::Sigma) && ((peak.sigma() <= 0.0)) )
    return false;


  const double xmin = is_step ? cont->lowerEnergy() : std::max( cont->lowerEnergy(), peak.mean() - 2.5*peak.sigma() );
  const double xmax = is_step ? cont->upperEnergy() : std::min( cont->upperEnergy(), peak.mean() + 2.5*peak.sigma() );

  shared_ptr<const SpecUtils::EnergyCalibration> cal = data->energy_calibration();
  assert( cal && cal->valid() );
  const shared_ptr<const vector<float>> &channel_energies = cal->channel_energies();
  assert( channel_energies && channel_energies->size() );
  if( !channel_energies )
    throw logic_error( "chi2_significance_test: invalid energy calibration" );

  const shared_ptr<const vector<float>> &gamma_counts = data->gamma_counts();
  assert( gamma_counts );

  const size_t num_spec_channel = cal->num_channels();
  const double lower_channelf = cal->channel_for_energy(xmin);
  const double upper_channelf = cal->channel_for_energy(xmax);

  const size_t begin_channel = std::min( static_cast<size_t>(std::max(lower_channelf, 0.0) ), num_spec_channel );
  const size_t end_channel = std::min( static_cast<size_t>(std::ceil(upper_channelf)), num_spec_channel );
  const size_t num_roi_channel = end_channel - begin_channel;

  assert( begin_channel < channel_energies->size() );
  assert( end_channel <= channel_energies->size() );

  // Check if above energy range, or if a really narrow ROI, if so, just accept it - whatever.
  if( (begin_channel >= channel_energies->size()) || (end_channel > channel_energies->size())
     || (end_channel < (begin_channel + 2)) )
  {
    return true;
  }

  const float * const energies = channel_energies->data() + begin_channel;
  const float * const channel_counts = gamma_counts->data() + begin_channel;

  vector<double> original_counts( num_roi_channel, 0.0 );

  for( const PeakDef &p : other_peaks )
    p.gauss_integral( energies, &(original_counts[0]), num_roi_channel );

  vector<double> without_peak_counts;

  // Build arrays of ROI peak pointers for offset_integral calls
  vector<const PeakDef *> all_roi_peak_ptrs;   // peak + other_peaks
  vector<const PeakDef *> other_peak_ptrs;     // just other_peaks (without the peak under test)

  all_roi_peak_ptrs.reserve( 1 + other_peaks.size() );
  all_roi_peak_ptrs.push_back( &peak );
  for( const PeakDef &p : other_peaks )
    all_roi_peak_ptrs.push_back( &p );

  other_peak_ptrs.reserve( other_peaks.size() );
  for( const PeakDef &p : other_peaks )
    other_peak_ptrs.push_back( &p );

// We can either fit the continuum with out the peak of interest, or just use the continuum as fit
//  I think the better thing to do is refit it, but this is totally untested - so the #define will
//  let you switch between options.
#define REFIT_CONTIUUM_FOR_NULL_HYPOTHESIS 1

  if( cont->type() == PeakContinuum::External )
  {
    shared_ptr<const Measurement> ext_contnuum = peak.continuum()->externalContinuum();
    for( size_t i = 0; !ext_contnuum && (i<other_peaks.size()); ++i )
      ext_contnuum = other_peaks[i].continuum()->externalContinuum();
    assert( ext_contnuum );
    if( !ext_contnuum || (ext_contnuum->num_gamma_channels() < data->num_gamma_channels()) )
      throw logic_error( "chi2_significance_test: invalid external continuum" );

    // In principle we should be able to just do a channel-by-channel sum, but we'll be safe
    //const float * const ext_counts = &(ext_contnuum->gamma_counts()->at(begin_channel));
    //for( size_t i = 0; i < num_roi_channel; ++i )
    //  original_counts[i] += ext_counts[i];

    cont->offset_integral( energies, &(original_counts[0]), num_roi_channel, data, all_roi_peak_ptrs.data(), all_roi_peak_ptrs.size() );
    without_peak_counts = original_counts;
  }//if( cont->type() == PeakContinuum::External )

  if( cont->type() != PeakContinuum::External )
  {
#if( REFIT_CONTIUUM_FOR_NULL_HYPOTHESIS )
    without_peak_counts = original_counts;
    const double ref_energy = cont->referenceEnergy();
    const double * const skew_pars = peak.coefficients() + static_cast<size_t>(PeakDef::CoefficientType::SkewPar0);

    // The null hypothesis is the same continuum family with the peak removed - not a different
    //  continuum model.  The peak-CDF step coefficients are bilinear with the peak amplitudes so
    //  the least-squares solve cannot fit them; hand it the fitted ones as known inputs, exactly
    //  as the fitters do.  Dropping the peak from the ROI already removes its own term from
    //  SUM_j(amp_j*CDFbar_j), which is the part of the step the peak is responsible for.
    const size_t num_poly_pars = PeakContinuum::num_linear_fit_pars( offset_type );
    const size_t num_step_pars = PeakContinuum::num_cdf_step_pars( offset_type );
    const vector<double> &fit_cont_pars = cont->parameters();
    assert( fit_cont_pars.size() == (num_poly_pars + num_step_pars) );
    const double * const step_coeffs = num_step_pars ? (fit_cont_pars.data() + num_poly_pars)
                                                     : nullptr;

    vector<double> amplitudes, continuum_coeffs, amp_uncerts, cont_uncerts;
    PeakFit::fit_amp_and_offset_imp(energies, channel_counts, nullptr, num_roi_channel, cont->type(),
                                    step_coeffs, ref_energy, {}, {}, other_peaks, peak.skewType(), skew_pars,
                                    amplitudes, continuum_coeffs, amp_uncerts, cont_uncerts, (double *)0 );

    // fit_amp_and_offset_imp returns only the polynomial coefficients; setParameters expects the
    //  step coefficients appended - and they are the fitted ones we just handed it as known.
    for( size_t k = 0; k < num_step_pars; ++k )
    {
      continuum_coeffs.push_back( fit_cont_pars[num_poly_pars + k] );
      cont_uncerts.push_back( 0.0 );
    }

    shared_ptr<PeakContinuum> tmp_continuum = make_shared<PeakContinuum>( *cont );
    tmp_continuum->setParameters( ref_energy, continuum_coeffs, cont_uncerts );
    tmp_continuum->offset_integral( energies, &(without_peak_counts[0]), num_roi_channel, data, other_peak_ptrs.data(), other_peak_ptrs.size() );

    cont->offset_integral( energies, &(original_counts[0]), num_roi_channel, data, all_roi_peak_ptrs.data(), all_roi_peak_ptrs.size() );
#else
    cont->offset_integral( energies, &(original_counts[0]), num_roi_channel, data, all_roi_peak_ptrs.data(), all_roi_peak_ptrs.size() );
    without_peak_counts = original_counts;
#endif
  }//if( cont->type() != PeakContinuum::External )


  peak.gauss_integral( energies, &(original_counts[0]), num_roi_channel );


  //Now compute Chi2
  double with_peak_chi2 = 0.0, without_peak_chi2 = 0.0;
  for( size_t i = 0; i < num_roi_channel; ++i )
  {
    const double ndata = channel_counts[i];
    const double uncert2 = (ndata > 1.0) ? ndata : 1.0;

    //TODO: This is ad-hoc, and just following PeakFitDiffCostFunction::parametersToPeaks implementation, for the moment
    if( ndata >= PEAK_FIT_MIN_CHANNEL_UNCERT )
    {
      const double diff_all_peaks = (ndata - original_counts[i]);
      with_peak_chi2 += (diff_all_peaks*diff_all_peaks) / uncert2;

      const double diff_without_peak = (ndata - without_peak_counts[i]);
      without_peak_chi2 += (diff_without_peak*diff_without_peak) / uncert2;
    }else
    {
      with_peak_chi2 += original_counts[i];
      without_peak_chi2 += without_peak_counts[i];
    }//if( ndata >= PEAK_FIT_MIN_CHANNEL_UNCERT ) / else
  }

  const double chi2Ratio = without_peak_chi2 / with_peak_chi2;

  const bool noDeltaRequired = (withoutPeakDSigma <= 0.0);
  bool noRatioRequired = (chi2ratioRequired <= 0.0);

  if( without_peak_chi2 < 5 )
    noRatioRequired = true;
  //Dont require the ratio test to apply if peaks share a continuum
  std::shared_ptr<const PeakContinuum> continuum = peak.continuum();
  for( const PeakDef &p : other_peaks )
    noRatioRequired |= (continuum == p.continuum());

  const double deltaChi2 = without_peak_chi2 - with_peak_chi2;


  static std::mutex s_mutex;
  //{
  //  std::lock_guard<std::mutex> lock( s_mutex );
  //  cerr << "chi2_significance_test(mean=" << peak.mean() << "): with_peak_chi2="
  //       << with_peak_chi2 << ", without_peak_chi2=" << without_peak_chi2 << endl
  //       << "chi2Ratio=" << chi2Ratio << " (requires " << chi2ratioRequired
  //       << "), deltaChi2=" << deltaChi2
  //       << " (required " << withoutPeakDSigma << ")"
  //       << ", noRatioRequired=" << noRatioRequired
  //       << ", peak.mean()=" << peak.mean()
  //       << ", peak.sigma()=" << peak.sigma()
  //       << " peak.amplitude()=" << peak.amplitude()
  //       << endl << endl;
  //}

  return ((noRatioRequired || (chi2Ratio >= chi2ratioRequired))
          && (noDeltaRequired || (deltaChi2 > withoutPeakDSigma)));
}//bool chi2_significance_test( ... 0

        
        









        
