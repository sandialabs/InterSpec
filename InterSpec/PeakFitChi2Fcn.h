#ifndef PeakFitChi2Fcn_h
#define PeakFitChi2Fcn_h
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
#include <memory>

#include "InterSpec/PeakDef.h"

#include "Math/IFunction.h"
#include "Minuit2/FCNBase.h"
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

namespace SpecUtils{ class Measurement; }


class PeakFitChi2Fcn
    : public ROOT::Minuit2::FCNBase
//      , public ROOT::Math::IBaseFunctionMultiDim //implemented, but not actually used, so left commented out
{
  //This class is a first attempt at better fitting than simple ROOT fitting,
  //  and is by no means correct or finished.
  //
  //ToDo:
  //  - The function chi2(...) should be gone through, refined, and checked that
  //    it is doing what it claims.
  //  - A -2*ln(liklihood) based figure of merrit method should be implemented.
  //  - Should fully integrate MultiPeakFitChi2Fcn into this class (along with
  //    checking of benchmarks)
  //  - Try inherriting from FCNGradientBase so can provide gradient

public:
  //FitPars gives the order of parameters.
  //  Currently extends the PeakDef::CoefficientType enum, in a kinda hacky way;
  //  there are some static asserts that try to ensure consistentcy incase
  //  PeakDef::CoefficientType is changed, hopefully to compile,
  //  PeakFitChi2Fcn::PeakPars will have to be altered as well.
  enum FitPars
  {
    Mean               = PeakDef::CoefficientType::Mean,
    Sigma              = PeakDef::CoefficientType::Sigma,
    GaussAmplitude     = PeakDef::CoefficientType::GaussAmplitude,
    SkewPar0           = PeakDef::CoefficientType::SkewPar0,
    SkewPar1           = PeakDef::CoefficientType::SkewPar1,
    SkewPar2           = PeakDef::CoefficientType::SkewPar2,
    SkewPar3           = PeakDef::CoefficientType::SkewPar3,
    Chi2DOF            = PeakDef::CoefficientType::Chi2DOF,
    ContinuumInfoField = PeakDef::CoefficientType::NumCoefficientTypes,
    SkewInfoField,
    RangeStartEnergy,
    RangeEndEnergy,
    OffsetEnergyRelativeTo,
    OffsetPolynomial0,
    OffsetPolynomial1,
    OffsetPolynomial2,
    OffsetPolynomial3,
    OffsetPolynomial4,
    NumFitPars
  };//enum PeakPars

public:
  //Constructor to fit full range of histograms
  PeakFitChi2Fcn( const int npeaks,
                  std::shared_ptr<const SpecUtils::Measurement> data,
                  std::shared_ptr<const SpecUtils::Measurement> continuum );

  //Constructor to fit partial range of histogram; only bins between lowerbin
  //  and upperbin will contribute to the chi2
  PeakFitChi2Fcn( const int npeaks,
                  const size_t lower_bin, const size_t upper_bin,
                  std::shared_ptr<const SpecUtils::Measurement> data,
                  std::shared_ptr<const SpecUtils::Measurement> continuum );

  void useReducedChi2( const bool use = true );        //default is true
  
  //setUseMultiPeakPunishment(bool): sets whether or not when fitting for
  //  multiple peaks that share the same continuum (via explicit shared
  //  polynomial continuum, or if sharing an External continuum), if the chi2
  //  should be increased if a peak isn't statistically significant, or
  //  if peaks are too near each other.
  //Not particularly well tested as of 20131230
  void setUseMultiPeakPunishment( bool punish );
  
  void estimateContinum( const bool estimate = true ); //default is false

  virtual unsigned int NDim() const;
  PeakFitChi2Fcn *Clone() const;
  PeakFitChi2Fcn &operator=( const PeakFitChi2Fcn &rhs );

  virtual double Up() const;

  virtual double operator()( const double *x ) const;
  virtual double operator()( const std::vector<double>& params ) const;  //does the work
  double chi2( const double *params ) const;

  double evalMultiPeakPunishment( const std::vector<const PeakDef *> &peaks ) const;
  
  //Creates the peaks from the given parameters, and if passed in the parameter
  //  errors, adds those to the peaks as well.
  //Note that the peaks will not have a Chi2DOF set, so you must do this
  //  manually.
  void parametersToPeaks( std::vector<PeakDef> &peaks,
                          const double *parameters,
                          const double *errors = 0 ) const;
  
  //addPeaksToFitter(...): adds peaks to the Minuit2::MnUserParameters,
  //  according to AddPeaksToFitterMethod specified.  This function should be
  //  used exclusively to add peaks to the MnUserParameters to ensure
  //  parametersToPeaks(...), and hence the chi2 function, will reconstruct the
  //  peaks correctly from the fit parameters.
  //This function will properly take into account peaks which share continuums,
  //  however, only for peaks passed in at one time.
  enum AddPeaksToFitterMethod
  {
    //All parameters fit for, assuming input is reasonable guess
    kFitForPeakParameters,
    
    //All parameters are held fixed
    kFixPeakParameters,
    
    //Mean, width, and skew parameters are held pretty close to current values
    kRefitPeakParameters,
    
    //Use when user indicated the location of the peak.  To indicate there is
    //  and external detector resolution defined, set the sigma error to > 0.0
    kFitUserIndicatedPeak
  };//enum AddPeaksToFitterMethod
  
  //addPeaksToFitter: throws exception if any of the input peaks are data
  //  defined
  static void addPeaksToFitter( ROOT::Minuit2::MnUserParameters &paramaters,
                                const std::vector<PeakDef> &near_peaks,
                                std::shared_ptr<const SpecUtils::Measurement> data,
                                AddPeaksToFitterMethod method );
  
  static void setSharedIndexToContinuumInfo( double &info, int index );
  static void setOffsetTypeToContinuumInfo( double &info,
                                            PeakContinuum::OffsetType type );
  
  static int continuumInfoToSharedIndex( double info );
  static PeakContinuum::OffsetType continuumInfoToOffsetType( double info );
  
  //testOffsetConversions(): tests that the continuum info successfully
  //  roundtrips to/from doubles (floating point math is tricky!)
//  static bool testOffsetConversions();
  
protected:
  bool m_useReducedChi2;
  bool m_useMultiPeakPunishment;
  size_t m_lower_channel, m_upper_channel;
  int m_npeaks;
  std::shared_ptr<const SpecUtils::Measurement> m_data;
  std::shared_ptr<const SpecUtils::Measurement> m_continium;
  
private:
  virtual double DoEval( const double *x ) const;
};//class PeakFitChi2Fcn



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
 ...X
 width     1               //m_numOffset + 0
 mean      1
 Amplitude 1  //If Amplitude==-999.9, then estimate from data height at peak centroid
 width     2  (multiple of first width, as fraction of total range; if zero, use first peak width; if negative use negative width)
 mean      2
 Amplitude 2  //If Amplitude==-999.9, then estimate from data height
 ...
 width     X  (if 0.0, then defer to coefficient at (m_numOffset+3); if negative, use negative value for width)
 mean      X
 Amplitude X  //If Amplitude==-999.9, then estimate from data height
 */

class MultiPeakFitChi2Fcn
: public ROOT::Minuit2::FCNBase
{
  public:
  static std::atomic<size_t> sm_ncalls;
  static std::atomic<bool> sm_call_opt_integrate;
  
  //This class is intended to fit for multiple peaks in a user defined region
  //  of the data, ignoring all other peaks.
  //  It currently does two things PeakFitChi2Fcn does not:
  //    1) Punishes for the shared peaks being to close
  //    2) Is optimized for speed of evaluation
public:
  /** Constructor to fit full range of histograms.
  */
  MultiPeakFitChi2Fcn( const int npeaks,
                      std::shared_ptr<const SpecUtils::Measurement> data,
                      PeakContinuum::OffsetType offsetType,
                      PeakDef::SkewType skewType,
                      const size_t lower_channel,
                      const size_t upper_channel );
  
  MultiPeakFitChi2Fcn &operator=( const MultiPeakFitChi2Fcn &rhs );
  
  
  virtual double Up() const;
  virtual double operator()( const std::vector<double>& params ) const;
  virtual double DoEval( const double *x ) const;
  double DoEval( const double *x, std::vector<PeakDef> &peaks ) const;
  
  void parametersToPeaks( std::vector<PeakDef> &peaks, const double *x,
                          const double *errors = 0 ) const;
  
  double dof() const;
  int nbin() const;
  
  void set_reldiff_punish_start( double reldiff_punish_start );
  void set_reldiff_punish_weight( double reldiff_punish_weight );
  
  //evalRelBinRange(...): returns chi2 of peaks for the data range.
  //  Note that bin numbers are internal bin numbers (so from 0 to nbin())
  double evalRelBinRange( const size_t beginRelChannel, const size_t endRelChannel,
                          const std::vector<PeakDef> &peaks ) const;
  
  //evalMultiPeakInsentive(...): returns punishment to the chi2 for peaks
  //  being to close, or for being statistically insignificant
  double evalMultiPeakInsentive( const std::vector<PeakDef> &peaks ) const;
  
protected:
  int m_npeak;
  size_t m_lower_channel, m_upper_channel;
  int m_numOffset;
  size_t m_nbin;
  
  double m_rangeLow, m_highRange;
  
  /** The lower channel energies for the channel range being used.
   Will be of size `m_nbin + 1` (i.e., the last entry is the upper energy of the last channel).
   The first entry in this vector is the lower energy of data channel `m_lower_channel`.
   */
  std::vector<float> m_energies;
  /** 
   
   Will be of size `m_nbin`.
   */
  std::vector<double> m_dataCounts;
  PeakContinuum::OffsetType m_offsetType;
  PeakDef::SkewType m_skewType;
  std::shared_ptr<const SpecUtils::Measurement> m_data;
  
  /* How close peaks have to be together before you start applying a weight
   to the chi2 to punish this - this keeps ROIs with more than one peak from
   collapsing to have multiple peaks located at the same mean and width.
   
   The Default value of 1.25 is good for situations where the user explicitly
   specifies how many peaks are in a ROI.
   A value of around 2.35482 is better when the number of peaks is being guessed
   in the fit.
   */
  double m_reldiff_punish_start;
  
  /* The weight to apply for punishment of near-by peaks.
   Default value of 2.0.
   A negative or zero value turns off this punishment.
   
   ToDo: Currently when this weight "turns on", the chi2 is very discontinous,
         need to define a much more reasonable weighting factor.
   */
  double m_reldiff_punish_weight;
  
  
  //If sm_call_opt_integrate is true, then will re-use the following peaks to
  //  avoid allocation overhead and such - currently only for development
  //  (looks to speed things up by maybe 8% or so?)
  mutable std::vector<PeakDef> m_workingpeaks;
};//class MultiPeakFitChi2Fcn



class LinearProblemSubSolveChi2Fcn
: public ROOT::Minuit2::FCNBase
{
  //This class uses a standard matrix based solution to solve for peak amplitude
  //  and continuum parameter values, so that minuit only has to solve for
  //  the means and width of the peaks.  In addition to this, this class also:
  //    1) Punishes for the peaks being to close
  //    2) Constrains the widths of the peaks to have a linearly related width
  //
  //For a single peak on a HPGE spectrum with a linear background, this chi2
  //  function takes about 45% as many minuit calls as PeakFitChi2Fcn, and
  //  produces similar results in a similar or less amount of cpu time.
  //
/*
   Want to fit for a number of peaks in the given range.
   All peaks should share a common polynomial continuum.
   Widths of the peaks are tied together.
   All peaks have same skew type, with same skew value(s).
   
   parameters:
   mean0     {mean of first peak}
   mean1     {mean of second peak}
   ...
   meanN     {mean of N'th peak}
   width0    {if npeaks==1, width of that peak; else width at m_lowerROI}
   widthFcn  {multiple of width0, plus fraction of total range time widthFcn}
   Skew0     {first skew parameter, if present - see `PeakDef::num_skew_parameters(SkewType)`}
   ....
   Skew3     {fourth, and last possible, skew parameter, if present}
*/

public:
  //Constructor most useful for when fitting new peaks (e.g. not refitting
  //  existing peaks).  The peaks created by the parametersToPeaks(...)
  //  function will not 'inherit' nuclide ID infor, or fit for values of any
  //  peviously existing peaks.
  LinearProblemSubSolveChi2Fcn( const size_t npeaks,
                               std::shared_ptr<const SpecUtils::Measurement> data,
                               const PeakContinuum::OffsetType offsetType,
                               const PeakDef::SkewType skewType,
                               const float lowerROI, const float upperROI );

  //Constructor useful when re-fitting peaks, to allow preserving of 'extra'
  //  information like nuclide/xray/reaction association, as well as accounting
  //  for cases when the input peaks specify to not fit for mean/sigma/amplitude
  //  of a peak when refitting.  It is assumed the order of peaks in
  //  'originalPeaks' is the same as the means in the params passed into
  //  operator()(...) or DoEval(...).  If the mean for one of the input peaks
  //  is specified to not be fit for, you must still specify a value for it in
  //  the parameters passed into operator()(...)/DoEval(...), although the
  //  value specified for in the input peak will be used.
  //  You should also becareful of the case where two peaks means could swap
  //  places in the fit, causing a potential confusion of nuclide ID (or
  //  whether to fit for values, etc), but I'm not convinced this is a huge
  //  problem/concern.
  LinearProblemSubSolveChi2Fcn(
        const std::vector< std::shared_ptr<const PeakDef> > &originalPeaks,
                               std::shared_ptr<const SpecUtils::Measurement> data,
                               const PeakContinuum::OffsetType offsetType,
                               const PeakDef::SkewType skewType,
                               const float lowerROI, const float upperROI );
  
  virtual double Up() const;
  virtual double operator()( const std::vector<double> &params ) const;
  virtual double DoEval( const double *x ) const;
  
  //parametersToPeaks(...): returns chi2 for the current paramters.  May throw
  //  if the linear sub problem cant be solved for.
  //If input peaks are specified when calling the LinearProblemSubSolveChi2Fcn
  //  contstructor, then output peaks will inherit nuclide associations, and
  //  fit for values.  If the input peaks specify to not fit for mean or width,
  //  this is enforced in the output peaks and calculation of chi2; if the
  //  amplitude is requested to not be fit for, this is reflected in the output
  //  peaks, but for chi2 computation a fit values of amplitude is used.
  //Output peaks will always be sorted accroding to mean.
  //The return value of this function is chi2, not chi2/dof.
  double parametersToPeaks( std::vector<PeakDef> &peaks,
                            const double *x,
                            const double *errors = 0 ) const;
  
  double dof() const;
  size_t nbin() const;
  size_t nfitPars() const;
  
  //punishment(..): the punishment to the chi2 for peaks being to close together
  //  or not statistically significant
  double punishment( const std::vector<PeakDef> &peaks ) const;
  double closenessPunishment( const std::vector<double> &means,
                              const std::vector<double> &sigmas ) const;
  
  /** Determines skew type, and starting parameters to use, form previous peaks in ROI.
   */
  static void skewFromPrevPeaks( const std::vector<std::shared_ptr<const PeakDef> > &prev_peaks,
                                 PeakDef::SkewType &skew_type,
                                 std::vector<bool> &fit_parameter,
                                 std::vector<double> &starting_value );
  
  /** You can use this function to add the skew fitting parameters to Minuit, when you dont have any previous
   values to use for the skew coefficients.
   */
  static void addSkewParameters( ROOT::Minuit2::MnUserParameters &pars,
                                 const PeakDef::SkewType skewType );
  
  static void addSkewParameters( ROOT::Minuit2::MnUserParameters &pars,
                                 const PeakDef::SkewType skewType,
                                const std::vector<bool> &fit_parameter,
                                const std::vector<double> &starting_value );
  
protected:
  void init( std::shared_ptr<const SpecUtils::Measurement> data );
  
protected:
  size_t m_nbin;
  const size_t m_npeak;
  const float m_lowerROI, m_upperROI;
  std::vector<float> m_y;
  std::vector<float> m_x;
  PeakContinuum::OffsetType m_offsetType;
  PeakDef::SkewType m_skewType;
  std::vector< std::shared_ptr<const PeakDef> > m_originalPeaks;
};//class LinearProblemSubSolveChi2Fcn



#endif
