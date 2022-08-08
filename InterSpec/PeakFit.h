#ifndef PeakFit_h
#define PeakFit_h
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

#include <deque>
#include <tuple>
#include <vector>

#include "Minuit2/FCNBase.h"

#include "InterSpec/PeakDef.h"


class DetectorPeakResponse;

namespace ROOT
{
  namespace Minuit2
  {
    class MnUserParameters;
  }//namespace Minuit2
}//namespace ROOT

namespace SpecUtils{ class Measurement; }


typedef std::shared_ptr<const DetectorPeakResponse> DetctorPtr;
typedef std::vector< std::shared_ptr<const PeakDef> > PeakShrdVec;

// TODO: put everything in this file into a namespace


struct SavitzyGolayCoeffs
{
  const int num_left;          //number of coef. to left of current point
  const int num_right;         //number of coef. to right of current point
  const int polynomial_order;  //order of smoothing polynomial,
                               //  (also highest conserved moment)
  const int ld;                //has to do with what derivative you want
  std::vector<double> coeffs;

  SavitzyGolayCoeffs( int bins_left, int bins_right,
                      int order, int derivative = 0 );
  void smooth( const std::vector<float> &input,
               std::vector<float> &output ) const;
  void smooth( const float *input, const int nSamples,
               std::vector<float> &output ) const;
};//struct SavitzyGolayCoeffs


//calculateContinuum(...) and findPeaksByRelaxation(...) are adapted
//  from ROOTS TSpectrum.  See implementation code for credits/details
enum
{ //An enum for calculateContinuum and findPeaksByRelaxation
  kBackOrder2 = 2, kBackOrder4 = 4, kBackOrder6 = 6, kBackOrder8 = 8,
  kBackIncreasingWindow, kBackDecreasingWindow,
  kBackSmoothing3, kBackSmoothing5, kBackSmoothing7,
  kBackSmoothing9, kBackSmoothing11, kBackSmoothing13,kBackSmoothing15
};//enum
const char *calculateContinuum( float *spectrum, int ssize,
                                       int numberIterations,
                                       int direction, int filterOrder,
                                       bool smoothing, int smoothWindow,
                                       bool compton );

//estimateContinuum(): estimates continuum for the data passed in using
//  "standard" parameters
std::shared_ptr<SpecUtils::Measurement> estimateContinuum( std::shared_ptr<const SpecUtils::Measurement> data );

std::vector<float> findPeaksByRelaxation( float *source, float *dest, int ssize,
                                          float sigma, double threshold,
                                          bool bckgrndRemove, int nIterations,
                                          bool markov, int averWindow );

std::vector< std::vector<PeakDef> >
         causilyDisconnectedPeaks( const double x0, const double x1,
                                   const double ncausality,
                                   const bool useRoiAsWell,
                                   const std::vector<PeakDef> &input_peaks );

std::vector< std::vector<std::shared_ptr<const PeakDef> > >
causilyDisconnectedPeaks(  const double ncausality,
                           const bool useRoiAsWell,
                           std::vector< std::shared_ptr<const PeakDef> > input_peaks );

//setPeakXLimitsFromData(...): intended to set how far to the right and left of
//  the mean a peaks continuum will be valid

/** Function to make sure the continuum of the defined peaks are unique to the peaks passed in.
 This function will
 
 This function is necassary when fitting new peaks that even though the PeakDef objects themselves
 get copied, so the original input peaks wont be modified, we need to make sure the continuum itself
 wont get modified as well.  The PeakDef tracks its continuum as a shared pointer, that may be
 shared by several PeakDefs, and if you copy a PeakDef, the pointer is just dumbly copied, meaning
 if you modify the continuum of a copied PeakDef, the continuum of the original PeakDef is also
 modified since they are the same object in memory.
 This is a poor design.  The continuum should own the PeakDef, not the other way around, but this
 function acts as a scab around this poor design, for the moment.

 */
void unique_copy_continuum( std::vector<PeakDef> &peaks );



//Note: smoothSpectrum(...) does not divide by bin widths
void smoothSpectrum( std::shared_ptr<const SpecUtils::Measurement> meas, const int side_bins,
                     const int order, const int derivative,
                     std::vector<float> &results );
void smoothSpectrum( const std::vector<float> &spectrum, const int side_bins,
                    const int order, const int derivative,
                    std::vector<float> &results );

#define PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL 0

//expected_peak_width_limits(): gives the extremes of expected peak
//  sigmas for detectors - the min sigma is 0.5 (0.75) times the smallest known
//  detector width, and the max sigma is 3 (2) times the largest known detector
//  width for highres (lowres) detectors.
//  Specify highres=false for NaI, CsI, LaBr3, or CZT detectors. lowres=true for
//  HPGe detectors.
void expected_peak_width_limits( const float energy,
                                 const bool highres,
                                 float &min_sigma_width_kev,
                                 float &max_sigma_width_kev );

//find_roi_for_2nd_deriv_candidate(): Takes the (smoothed) second derivative
//  and finds the second intersection of it with y==0, on each side of
//  'peakmean', then for this region, it finds the line which just touches the
//  (smoothed) data, such that all channel content heights in between is above
//  this line.  The two points where this line touches the data define the
//  initial estimates for 'lowerEnengy' and 'upperEnergy'.  The extent beyond
//  the mean will try to be extended on each side by a maximum of 35%, if the
//  line stays within 5 sigma (cumulative) of the data - this is to agree a
//  little better with what an analysts eyes would choose, and isnt super
//  necassary.
//Throws exception on error, for instance if the second derivative at 'peakmean'
//  isnt negative, or data is invalid.
void find_roi_for_2nd_deriv_candidate( double &lowerEnengy, double &upperEnergy,
                            const float peakmean,
                            const std::shared_ptr<const SpecUtils::Measurement> &data );

//For meaning of stat_threshold and hypothesis_threshold see notes for
//  the chi2_significance_test(..) function
//The fixedpeaks passed in are not included in the results, and are taken as
//  fixed in the fit, and will not be deleted even if they are not significant.
//If isRefit==true, then tighter restrictions are made on the mean and widths
//  of the peaks; see kRefitPeakParameters notes above.
//Results includes all the 'all_peaks' passed in, with the ones in the specified
//  range having been refit.
std::vector<PeakDef> fitPeaksInRange( const double x0, const double x1,
                                      const double ncausalitysigma,
                                      const double stat_threshold,
                                      const double hypothesis_threshold,
                                      std::vector<PeakDef> all_peaks,
                                      std::shared_ptr<const SpecUtils::Measurement> data,
                                      const std::vector<PeakDef> &fixedpeaks,
                                      bool isRefit = false );



//secondDerivativePeakCanidatesWithROI(): similar to above function,
//
//  start_bin, end_bin: allow you to specify a the search range; if a smaller
//  or equal channel number is specified for end_channel than start_channel,
//  or a to large of number is specified, then dataH->num_gamma_channels() will
//  be used.
std::vector<std::shared_ptr<PeakDef> > secondDerivativePeakCanidatesWithROI(
                                                          std::shared_ptr<const SpecUtils::Measurement> data,
                                                           size_t start_channel,
                                                           size_t end_channel );

/** Similar to #secondDerivativePeakCanidatesWithROI, but doesnt waste time finding
 the ROI for each peak.  Provides reults as a tuple of {mean,sigma,area} for
 each found peak.
 
 Takes about 40% of the time as #secondDerivativePeakCanidatesWithROI
 */
void secondDerivativePeakCanidates( const std::shared_ptr<const SpecUtils::Measurement> data,
                                    size_t start_channel,
                                    size_t end_channel,
                                   std::vector< std::tuple<float,float,float> > &results );


//combine_peaks_to_roi: throws exception on error
void combine_peaks_to_roi( PeakShrdVec &coFitPeaks,
                          double &roiLower,
                          double &roiUpper,
                          bool &lowstatregion,
                          const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                          const PeakShrdVec &inpeaks,
                          const double mean0,
                          const double sigma0,
                          const double area0,
                          const double pixelPerKev );

bool check_lowres_single_peak_fit( const std::shared_ptr<const PeakDef> peak,
                                  const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                  const bool lowstatregion,
                                  const bool automated );

void get_candidate_peak_estimates_for_user_click(
                                                 double &sigma0, double &mean0, double &area0,
                                                 const double x,
                                                 const double pixelPerKev,
                                                 const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                                 const PeakShrdVec &inpeaks );

//searchForPeakFromUser: looks for a peak near x.  The returned pair contains
//  peaks that should be added (in .first) and peaks that should be removed
//  (in .second)
//Tries to take into account user resolution in limiting how far the fit peak
//  can be away from the nominal energy.
//  If pixelPerKev <= 0.0 is specified, then it is assumed this is an automated
//  search, and tougher quality requirements will be placed on the fit peaks.
std::pair< PeakShrdVec, PeakShrdVec > searchForPeakFromUser( const double x,
                                           double pixelPerKev,
                                           const std::shared_ptr<const SpecUtils::Measurement> &data,
                                           const PeakShrdVec &existing_peaks,
                                           std::shared_ptr<const DetectorPeakResponse> drf );

//refitPeaksThatShareROI: intended to refit peaks fit for by
//  searchForPeakFromUser(...), for instance when you modfiy the ROI range.
//  Returns an empty result if error occurs, or is not able to produce a better
//  fit
//  meanSigmaVary is how many sigma we should limit the peak means to. A negative
//  value means no limit (untested).  Recomend 0.25.
std::vector< std::shared_ptr<const PeakDef> >
    refitPeaksThatShareROI( const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                            const DetctorPtr &detector,
                            const std::vector< std::shared_ptr<const PeakDef> > &inpeaks,
                            const double meanSigmaVary );


//For meaning of stat_threshold and hypothesis_threshold see notes for
//  the chi2_significance_test(..) function
//The fixedpeaks passed in are not included in the results, and are taken as
//  fixed in the fit, and will not be deleted even if are not significant;
//  also, they do not influence the X-range fit for.
void fitPeaks( const std::vector<PeakDef> &input_peaks,
                      const double stat_threshold,
                      const double hypothesis_threshold,
                      std::shared_ptr<const SpecUtils::Measurement> data,
                      std::vector<PeakDef> &results,
                      const std::vector<PeakDef> &fixedpeaks,
                      bool amplitudeOnly ) throw();

enum MultiPeakInitialGuesMethod
{
  UniformInitialGuess,    //With a rough prefit performed
  FromDataInitialGuess,   //Tries to use second derivative
  FromInputPeaks          //uses peaks populated in 'answer', throws if answer.size!=nPeaks, or if they dont all share a continuum
};//enum MultiPeakInitialGuesMethod

//findPeaksInUserRange(...): method to simultaneously fit for 'nPeaks' in the
//  the range x0 to x1
void findPeaksInUserRange( double x0, double x1, int nPeaks,
                          MultiPeakInitialGuesMethod method,
                          std::shared_ptr<const SpecUtils::Measurement> dataH,
                          std::shared_ptr<const DetectorPeakResponse> detector,
                          std::vector<std::shared_ptr<PeakDef> > &answer,
                          double &chi2 );


void findPeaksInUserRange_linsubsolve( double x0, double x1, int nPeaks,
                          MultiPeakInitialGuesMethod method,
                          std::shared_ptr<const SpecUtils::Measurement> dataH,
                          std::shared_ptr<const DetectorPeakResponse> detector,
                          std::vector<std::shared_ptr<PeakDef> > &answer,
                          double &chi2 );

//nsigma is number of sigma away from mean of gaussian being tested to
//  consider, inorder for the peak to be in the range; for non-gaus peaks
//  peak.lowerX() and peak.upperX() are used regardless of nsigma.
//  returned peaks are sorted by mean.
std::vector<PeakDef> peaksInRange( const double lowx,
                                   const double highx,
                                   const double nsigma,
                                   const std::vector<PeakDef> &inputs );

//peaksTouchingRange(...): similar to peaksInRange(...), but instead uses the
//  peaks own definition of the peaks range (e.g. lowerX(), upperX()) to decide
//  if the peaks are in the specified range.
std::vector<PeakDef> peaksTouchingRange( double lowx, double highx,
                                         const std::vector<PeakDef> &inputs );
PeakShrdVec peaksTouchingRange( double lowx, double highx,
                                const PeakShrdVec &inputs );


//fit_to_polynomial(...):  Fits data to a polynomial function using Linear Least
//  Squares, returning the chi2 of the fit.  Assumes data is Poisson distributed.
//  'polynomial_order' 0 is a constant, 1 is a straight line, and so on.
//  Implemenetation is not very computationally efficient.
//  To calculate the y for a given x, using the result you would:
//    double y = 0.0;
//    for( int i = 0; i <= polynomial_order; ++i )
//      y += poly_coeffs[i] * std::pow( x, double(i) );
//  or call evaluate_polynomial(...)
//Will throw exception if run into numerical issues
double fit_to_polynomial( const float *x, const float *y, const size_t nbin,
                         const int polynomial_order,
                         std::vector<double> &poly_coeffs,
                         std::vector<double> &coeff_uncerts );

//evaluate_polynomial(): evaluates the polynomia found by fit_to_polynomial(...)
//  for the given x value and coefficients
double evaluate_polynomial( const double x,
                           const std::vector<double> &poly_coeffs );


/** Fits the continuum and amplitude of peaks with specified means and sigmas, over the data range specified.  Uses a
 matrix based linear regression fitter to perform minization.
 
 @param energies The lower-channel energies of ROI.  ROI defined by energies[0] to energies[nbin]. Must be of at least length nbin+1.
 @param data The channel counts of the ROI.  Must be of at least length nbin.
 @param nbin The number of channels in the ROI.
 @param num_polynomial_terms The number of polynomial continuum terms to fit for.
        0 is no continuum (untested), 1 is constant, 2 is linear sloped continuum, etc
 @param step_continuum Specifies whether or not a step in the continuum is to be used.  Currently
 
 Currently the implementation is reasonably inefficient.
 
 Throws exception upon ill-posed input.
 */
double fit_amp_and_offset( const float *energies, const float *data, const size_t nbin,
                           const int num_polynomial_terms,
                           const bool step_continuum,
                           const double ref_energy,
                           const std::vector<double> &means,
                           const std::vector<double> &sigmas,
                           const std::vector<PeakDef> &fixedAmpPeaks,
                           std::vector<double> &amplitudes,
                           std::vector<double> &continuum_coeffs,
                           std::vector<double> &amplitudes_uncerts,
                           std::vector<double> &continuum_coeffs_uncerts );

/** Get the chi2 and degrees of freedom for a peaks that share a ROI.
 All peaks in 'peaks' _must_ share the same continuum (or else assert will happen).
 */
void get_chi2_and_dof_for_roi( double &chi2, double &dof,
                               const std::shared_ptr<const SpecUtils::Measurement> &data,
                               const std::vector<PeakDef *> &peaks );

//set_chi2_dof(): computes and sets the Chi2/Dof for gaussian peaks with index
// 'startpeakindex' through 'startpeakindex + npeaks'.
//  Takes into account sharing of ROI between peaks.
//  returns the total number of degrees of freed, of all ROIs
double set_chi2_dof( std::shared_ptr<const SpecUtils::Measurement> data,
                   std::vector<PeakDef> &fitpeaks,
                   const size_t startpeakindex, const size_t npeaks );

//chi2_for_region(...): gives the chi2 or a region of data, given
//  the input peaks.
double chi2_for_region( const PeakShrdVec &peaks,
                        const std::shared_ptr<const SpecUtils::Measurement> &data,
                        const int lowBin,
                        const int highBin );

//stat_threshold: this is how incompatible with background/continuum the data
//                must be, before a peak is allowed to exist.  Reasonable
//                numbers for this are probably between 1 and 5.
//hypothesis_threshold: this specifies how well the peak must match in shape
//                      to a gaussian in order to keep the peak.  The higher
//                      this number, the more like a gaussian it must be. It
//                      is actually the ratio of the null hypothesis chi2
//                      to the test hypothesis chi2.  A reasonable value for
//                      this seems to be 4.
bool chi2_significance_test( PeakDef peak,
                             const double stat_threshold,  //this is how large the chi2 without the peak must be, inorder for peak to be necessary - the higher the number this is, the further away from the background the data must be before a peak becomes necaassary
                             const double hypothesis_threshold,  //this says roughly how good a gaussian explains the excess of data over background - higher the number input, the more closer to a gaussian the shape has to be
                             std::vector<PeakDef> other_peaks,
                             std::shared_ptr<const SpecUtils::Measurement> data );


namespace ExperimentalAutomatedPeakSearch
{
  std::vector<std::shared_ptr<const PeakDef> >
              search_for_peaks( const std::shared_ptr<const SpecUtils::Measurement> meas,
                                const std::shared_ptr<const DetectorPeakResponse> drf,
                                std::shared_ptr<const std::deque< std::shared_ptr<const PeakDef> > > origpeaks,
                                const bool singleThreaded );
}//namespace ExperimentalAutomatedPeakSearch


#define WRITE_CANDIDATE_PEAK_INFO_TO_FILE 0
#define WRITE_CANDIDATE_PEAK_TERMINAL_DEBUG_LEVEL 0

namespace ExperimentalPeakSearch
{
  
  //search_for_peaks(): a convienience function to call below
  //  search_for_peaks(...) that uses the current best guess of arguments
  std::vector<PeakDef> search_for_peaks( const std::shared_ptr<const SpecUtils::Measurement> meas,
                                         const std::vector<PeakDef> &origpeaks );
  
  
bool find_spectroscopic_extent( std::shared_ptr<const SpecUtils::Measurement> meas,
                               size_t &lower_channel,
                               size_t &upper_channel );

std::vector<PeakDef> search_for_peaks( const std::shared_ptr<const SpecUtils::Measurement> meas,
                                      const double min_chi2_dof_thresh,
                                      const double min_gross_counts_sig_thresh,
                                      const double above_line_chi2_thresh,
                                      const int side_bins,
                                      const int smooth_order,
                                      const double second_deriv_thresh,
                                      const double stat_thresh,
                                      const double width_thresh,
                                      const std::vector<PeakDef> &origpeaks /*included in result, unmodified, wont have duplciate */
#if( WRITE_CANDIDATE_PEAK_INFO_TO_FILE )
                                      , std::shared_ptr<const DetectorPeakResponse> detector
#endif
                                       );
  
  
  class AutoPeakSearchChi2Fcn : public ROOT::Minuit2::FCNBase
  {
  public:
    enum ResolutionType
    {
      Polynomial0thOrder,
      Polynomial1stOrder,
      Polynomial2ndOrder,
      Polynomial3rdOrder,
      SqrtEnergy
    };//enum ResolutionType
    
  protected:
    bool m_inited;
    std::shared_ptr< const std::vector<float> > m_y;
    std::shared_ptr< const std::vector<float> > m_x;
    std::shared_ptr<const SpecUtils::Measurement> m_meas;
    
    std::vector<PeakDef> m_candidates;
    std::vector<PeakDef> m_fixed_peaks;
    
    std::vector< std::vector<bool> > m_grouped_isfixed;
    std::vector< std::vector<PeakDef> > m_grouped_candidates;
    
    //precalculate the begining and ending iterators over energy and counts
    std::vector< std::pair<const float *,const float *> > m_group_counts;
    std::vector< std::pair<const float *,const float *> > m_group_energies;
    
    size_t m_lower_channel;
    size_t m_upper_channel;
    
    size_t m_num_peakstotal;
    size_t m_nbins_used;
    
    ResolutionType m_resolution_type;
    
  public:
    int m_side_bins;
    int m_smooth_order;
    double m_second_deriv_thresh;
    double m_stat_thresh;
    double m_width_thresh;
    
    //m_nsigma_near_group: how near, in terms of average sigma, two peaks must be
    //  to eachother in order to be in the same ROI.  The second-derivative ROI
    //  must also overlap as well.
    double m_nsigma_near_group;
    
    //m_min_chi2_dof_thresh: the minimum a peak must reduce a regions chi2/DOF
    //  in order to be kept by significance_test(...)
    double m_min_chi2_dof_thresh;
    double m_min_gross_counts_sig_thresh;
    
  public:
    //fixed_peaks: peaks that will be subtracted from data before finding aditional
    //             candidate peaks.  These peaks will then be included in the fit,
    //             and their properties allowed to vary, but in a more restricted
    //             manor then the new peaks.
    AutoPeakSearchChi2Fcn( std::shared_ptr<const SpecUtils::Measurement> data,
                          const std::vector<PeakDef > &fixed_peaks );    
    
    size_t lower_spectrum_channel() const;
    size_t upper_spectrum_channel() const;
    size_t num_initial_candidate_peaks() const;
    
    void second_derivative( const std::vector<float> &input, std::vector<float> &results );
    
    std::vector<PeakDef> candidate_peaks( const std::vector<float> &energies,
                                         const std::vector<float> &channel_counts );
    
    
    
    //init(): returns if fitting can proced;
    bool init();
    
    virtual double Up() const;
    virtual double operator()( const std::vector<double> &params ) const;
    double peak_sigma( const double energy, const std::vector<double> &pars ) const;
    
    ROOT::Minuit2::MnUserParameters initial_parameters() const;
    
    double eval_chi2( const std::vector<double> &params ) const;
    
    void fit_peak_group( const std::vector<PeakDef> &peaks,
                        const std::vector<double> &resolution_coefs,
                        const size_t group,
                        const double *pars,
                        double &chi2,
                        std::vector<PeakDef> &results ) const;
    //pars_to_peaks(): returns chi2
    double pars_to_peaks( std::vector<PeakDef> &resultpeaks,
                         const std::vector<double> &pars,
                         const std::vector<double> &errors = std::vector<double>() ) const;
    
    //punishment(..): the punishment to the chi2 for peaks being to close together
    //  or not statistically significant
    double punishment( const std::vector<PeakDef> &peaks ) const;
    
    double closeness_punishment( const std::vector<double> &means,
                                const std::vector<double> &sigmas ) const;
    //other_peaks may contain peak
    bool significance_test( const PeakDef &peak,
                           std::vector<PeakDef> other_peaks,
                           double *chi2DOF = 0,
                           double *grosCountsSigma = 0
                           ) const;
    
    std::vector<PeakDef> filter_peaks( const std::vector<PeakDef> &all_peaks ) const;
  };//class AutoPeakSearchChi2Fcn

  

}//namespace ExperimentalPeakSearch

#endif //#ifndef PeakFit_h
