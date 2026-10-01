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


#ifdef _MSC_VER
#undef isinf
#undef isnan
#undef isfinite
#undef isnormal
#endif

#include "ceres/ceres.h"

/** How a ROI is fit by Poisson maximum likelihood - by default only the sparse ones (see
 `sm_sparse_data_likelihood_threshold`), or every ROI with `PeakFitLMOptions::ForcePoissonLikelihood`:
   0: IRLS - after the chi2 fit, its solve is repeated with each channel's variance frozen at the
      previous model, with a deviance line search, until the peaks stop changing.
   1: all-in-Ceres - amplitudes and continuum become Ceres parameters, and the Poisson deviance is
      minimized starting from the chi2 fit.  Only applies through `fit_peaks_in_roi_LM(...)`, and not
      to an External continuum; a joint fit of several ROIs (`fit_peaks_in_spectrum_LM(...)` with a
      skew type) still uses IRLS.
 The two give the same areas and calibrated uncertainties for peaks of significance z >= 2 (2026-10
 evaluation); all-in-Ceres takes ~1.5x chi2's CPU against IRLS's 2-4x, while IRLS has the smaller
 tail of poor fits, is better for the weakest peaks, and stays closer to chi2 when the peak model
 does not describe the data.
 */
#define SPARSE_DATA_LIKELIHOOD_USE_CERES 0

#include "InterSpec/PeakDef.h"
#include "InterSpec/PeakFit.h"
#include "SpecUtils/SpecFile.h"
#include "InterSpec/PeakDists.h"
#include "InterSpec/PeakFitLM.h"

#if( PEAK_FIT_LM_PARALLEL_ROIS )
#include <future>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include "SpecUtils/SpecUtilsAsync.h"
#endif
#include "InterSpec/PeakFitUtils.h"
#include "InterSpec/PeakFitDetPrefs.h"
#include "InterSpec/PeakFitSpecImp.h"
#include "SpecUtils/EnergyCalibration.h"
#include "InterSpec/DetectorPeakResponse.h"

#include "InterSpec/PeakFit_imp.hpp"
#include "InterSpec/PeakDists_imp.hpp"
#include "InterSpec/RelActCalcAuto_imp.hpp"
#include "InterSpec/PeakFitLMObjective_imp.hpp"

// Undefine isnan and isinf macros
#undef isnan
#undef isinf


using namespace std;

namespace
{
/** Ceres iteration callback that aborts the solve when the cancellation token
 flips to true.  Holds a `shared_ptr` to the atomic so the flag is kept alive
 even if the originating `SpecMeas` is destroyed mid-fit.
 */
class CancelIterationCallback : public ceres::IterationCallback
{
  std::shared_ptr<const std::atomic<bool>> m_cancel_flag;

public:
  explicit CancelIterationCallback( std::shared_ptr<const std::atomic<bool>> flag )
    : m_cancel_flag( std::move(flag) )
  {
  }

  ceres::CallbackReturnType operator()( const ceres::IterationSummary & ) override
  {
    if( m_cancel_flag && m_cancel_flag->load( std::memory_order_relaxed ) )
      return ceres::SOLVER_ABORT;
    return ceres::SOLVER_CONTINUE;
  }
};//class CancelIterationCallback


/** To make the code of `PeakFitDiffCostFunction` work with either `PeakDef`
 or `RelActCalcAuto::PeakDefImp<T>`, we'll add in a few interface functions
 */
template<typename T>
struct PeakContinuumImp : public RelActCalcAuto::PeakContinuumImp<T>
{
  std::shared_ptr<const SpecUtils::Measurement> m_external_continuum;
  std::shared_ptr<const SpecUtils::Measurement> externalContinuum() const { return m_external_continuum; }
  void setExternalContinuum( const std::shared_ptr<const SpecUtils::Measurement> &data ){ m_external_continuum = data; }
};

template<typename T>
struct PeakDefImpWithCont : public RelActCalcAuto::PeakDefImp<T>
{
  std::shared_ptr<PeakContinuumImp<T>> m_continuum;
  std::shared_ptr<PeakContinuumImp<T>> getContinuum(){ return m_continuum; }
  void setContinuum( const std::shared_ptr<PeakContinuumImp<T>> &continuum ) { m_continuum = continuum;}
  void setMeanUncert( const T &mean_uncert ){ }
  void setAmplitudeUncert( const T &amp_uncert ){ }
  void setSigmaUncert( const T &sigma_uncert ){ }
  void inheritUserSelectedOptions( const PeakDef &parent, const bool inheritNonFitForValues [[maybe_unused]] )
  {
    //We dont really need this function, since we never use any of this info, but whatever.
    this->m_parent_nuclide = parent.parentNuclide();
    this->m_transition = parent.nuclearTransition();
    this->m_rad_particle_index = parent.decayParticleIndex();
    this->m_xray_element = parent.xrayElement();
    this->m_reaction = parent.reaction();
    this->m_src_energy = (this->m_parent_nuclide || this->m_xray_element || this->m_reaction)
                          ? parent.gammaParticleEnergy() : 0.0f;
    this->m_gamma_type = parent.sourceGammaType();
    //m_rel_eff_index = ...
  }

  static bool lessThanByMean( const PeakDefImpWithCont<T> &lhs, const PeakDefImpWithCont<T> &rhs )
  {
    return (lhs.m_mean < rhs.m_mean);
  }
};//struct PeakDefImpWithCont


/** Replaces contents of input peaks vector with completely new peaks, that crucually have had new continuums allocated. */
void local_unique_copy_continuum( vector<shared_ptr<const PeakDef>> &input_peaks )
{
  // Use group_peaks_by_roi for deterministic ordering (sorted by lowerEnergy)
  std::vector<std::pair<const PeakContinuum *, std::vector<std::shared_ptr<const PeakDef>>>> groups
    = group_peaks_by_roi( input_peaks );

  input_peaks.clear();
  for( auto &group : groups )
  {
    std::shared_ptr<PeakDef> first_copy = make_shared<PeakDef>( *group.second[0] );
    first_copy->makeUniqueNewContinuum();
    std::shared_ptr<PeakContinuum> newcont = first_copy->continuum();
    input_peaks.push_back( first_copy );

    for( size_t i = 1; i < group.second.size(); ++i )
    {
      std::shared_ptr<PeakDef> copy = make_shared<PeakDef>( *group.second[i] );
      copy->setContinuum( newcont );
      input_peaks.push_back( copy );
    }
  }
  std::sort( begin(input_peaks), end(input_peaks), &PeakDef::lessThanByMeanShrdPtr );
}//unique_copy_continuum(...)
}//namespace


namespace PeakFitLM
{

/** The statistic `PeakFitDiffCostFunction` minimizes.  IRLS is not an objective of its own: it is the
 chi2 objective with per-channel variances set between solves (see `run_ceres_fit(...)`).
 */
enum class FitObjective : int
{
  NeymanChi2,

  /** All-in-Ceres Poisson maximum likelihood (`SPARSE_DATA_LIKELIHOOD_USE_CERES`): amplitudes and
   continuum are Ceres parameters, and the residuals are the signed-root Poisson deviance. */
  PoissonAllInCeres
};//enum class FitObjective


static thread_local FitObjectiveDiagnostics tl_fit_objective_diagnostics;

FitObjectiveDiagnostics take_fit_objective_diagnostics()
{
  const FitObjectiveDiagnostics answer = tl_fit_objective_diagnostics;
  tl_fit_objective_diagnostics = FitObjectiveDiagnostics();
  return answer;
}


/** The options to actually fit with: when the data has a negative channel (e.g., a
 background-subtracted spectrum), where Poisson statistics do not apply, the likelihood is turned off.
 Idempotent.  Throws if both `ForcePoissonLikelihood` and `NoSparseDataLikelihood` are given.
 */
static Wt::WFlags<PeakFitLMOptions> effective_options( Wt::WFlags<PeakFitLMOptions> options,
                                                       const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  const bool forced = options.test( PeakFitLMOptions::ForcePoissonLikelihood );
  if( options.test( PeakFitLMOptions::NoSparseDataLikelihood ) )
  {
    if( forced )
      throw std::logic_error( "PeakFitLM: ForcePoissonLikelihood and NoSparseDataLikelihood are exclusive." );
    return options;
  }

  const std::shared_ptr<const std::vector<float>> counts = data ? data->gamma_counts() : nullptr;
  const bool has_negative = counts && std::any_of( begin(*counts), end(*counts), []( const float c ){ return c < 0.0f; } );
  if( !has_negative )
    return options;

  options.clear( PeakFitLMOptions::ForcePoissonLikelihood );
  options |= PeakFitLMOptions::NoSparseDataLikelihood;
  if( forced )
    tl_fit_objective_diagnostics.fell_back_to_chi2 = true;

  return options;
}//effective_options(...)


/** `options` for a plain chi2 fit: the likelihood off. */
static Wt::WFlags<PeakFitLMOptions> chi2_only_options( Wt::WFlags<PeakFitLMOptions> options )
{
  options.clear( PeakFitLMOptions::ForcePoissonLikelihood );
  return options | PeakFitLMOptions::NoSparseDataLikelihood;
}


/** For the Poisson-deviance residuals: expected counts below this fraction of the ROI's mean counts
 per channel (or of one count, if larger) are smoothly floored; as `min_expected_channel_counts(...)`
 in DetectionLimitCalc.
 */
const double sm_deviance_min_expected_frac = 1.0E-8;


/** A ROI is refit by Poisson maximum likelihood, by default, when
   sqrt( SUM 1/max(model,1) )
 over the channels of its peak region (within 1.5 FWHM of a peak mean) exceeds this; see
 `sparse_data_statistic(...)`.  The modified-Neyman chi2 is biased by about one count per channel, so
 for a level fit over N channels of mu counts the bias, in units of its uncertainty, is about
 sqrt(N/mu) - which is what the statistic estimates, channel by channel.

 Chosen with target/peak_fit_improve_ai's peak_fit_objective_eval (2026-10):
  - synthetic spectra with exact truth (26650 paired fits): below 1.5 the chi2 and maximum-likelihood
    areas differ systematically by less than 0.1 sigma and both are unbiased; above it chi2's bias
    grows to over 1 sigma while the likelihood fit stays unbiased.
  - GADRAS-inject HPGe spectra (~46k fits), against truth: the likelihood fit's median area error is
    the smaller from about 1 up, but below ~2.5 it also follows the (mild) peak-shape mismatch
    further than chi2 does, so its mean bias there is a little worse; above 2.5 it is better on both.
 2.0 takes the likelihood fit where both data sets favour it; about a third of the synthetic grid
 and half the inject fits (both deliberately low-statistics heavy) are above it.
 */
const double sm_sparse_data_likelihood_threshold = 2.0;


/** The sparse-data statistic described at `sm_sparse_data_likelihood_threshold`, for one ROI.
 `energies` holds `nchannel + 1` channel edges, `model` the fitted counts of each channel, and
 `mean_fwhms` the (mean, FWHM) of the ROI's peaks.
 */
static double sparse_data_statistic( const float * const energies, const double * const model, const size_t nchannel,
                                     const std::vector<std::pair<double,double>> &mean_fwhms )
{
  double sum_inverse = 0.0;
  for( size_t ch = 0; ch < nchannel; ++ch )
  {
    const double energy = 0.5*(static_cast<double>(energies[ch]) + static_cast<double>(energies[ch+1]));
    bool in_peak_region = false;
    for( const std::pair<double,double> &mf : mean_fwhms )
      in_peak_region |= (fabs( energy - mf.first ) <= 1.5*mf.second);
    if( in_peak_region )
      sum_inverse += 1.0 / std::max( model[ch], 1.0 );
  }
  return std::sqrt( sum_inverse );
}//sparse_data_statistic(...)


/** The model counts (continuum + peaks) of the `nchannel` channels starting at `ch0`, for peaks that
 share one continuum.
 */
static std::vector<double> roi_model_counts( const std::vector<std::shared_ptr<const PeakDef>> &peaks,
                                             const std::shared_ptr<const SpecUtils::Measurement> &data,
                                             const size_t ch0, const size_t nchannel )
{
  const float * const energies = data->channel_energies()->data() + ch0;
  std::vector<const PeakDef *> roi_peaks;
  for( const std::shared_ptr<const PeakDef> &p : peaks )
    roi_peaks.push_back( p.get() );

  std::vector<double> model( nchannel, 0.0 );
  peaks.front()->continuum()->offset_integral( energies, model.data(), nchannel, data,
                                               roi_peaks.data(), roi_peaks.size() );
  for( const std::shared_ptr<const PeakDef> &p : peaks )
    p->gauss_integral( energies, model.data(), nchannel );

  return model;
}//roi_model_counts(...)


/** `sparse_data_statistic(...)` for fitted peaks that share one continuum, from their model. */
static double sparse_data_statistic( const std::vector<std::shared_ptr<const PeakDef>> &peaks,
                                     const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  if( peaks.empty() || !data || !data->channel_energies() )
    return 0.0;

  const std::shared_ptr<const PeakContinuum> cont = peaks.front()->continuum();
  const size_t ch0 = data->find_gamma_channel( static_cast<float>( cont->lowerEnergy() ) );
  const size_t ch1 = data->find_gamma_channel( static_cast<float>( cont->upperEnergy() ) );
  if( ch1 < ch0 )
    return 0.0;

  const size_t nchannel = ch1 - ch0 + 1;
  const std::vector<double> model = roi_model_counts( peaks, data, ch0, nchannel );
  std::vector<std::pair<double,double>> mean_fwhms;
  for( const std::shared_ptr<const PeakDef> &p : peaks )
    mean_fwhms.emplace_back( p->mean(), p->fwhm() );

  return sparse_data_statistic( data->channel_energies()->data() + ch0, model.data(), nchannel, mean_fwhms );
}//sparse_data_statistic( peaks )


/** Smallest per-channel variance the reweighted (IRLS) objective uses, so a channel whose model is
 essentially zero cannot dominate the fit; relative to the ROI's mean counts per channel.
 */
const double sm_irls_min_variance_frac = 1.0E-3;

/** Maximum number of reweighted (IRLS) passes after the initial chi2 fit. */
const size_t sm_max_irls_passes = 12;

/** A reweighted pass that lowers the Poisson deviance by less than this has converged: a deviance
 change of 1E-4 corresponds to moving the parameters by about 0.01 of their uncertainty.
 */
const double sm_irls_deviance_tolerance = 1.0E-4;

/** A likelihood refit of a ROI already at its likelihood minimum gets back there by way of the chi2
 solution, so its deviance can come back above the input's by about `sm_irls_deviance_tolerance`;
 `refitPeaksThatShareROI_LM(...)` only counts a refit as worse than its input beyond this.
 */
const double sm_refit_deviance_tolerance = 10.0*sm_irls_deviance_tolerance;


/** Info about a single ROI (region of interest), used internally by PeakFitDiffCostFunction.
 A ROI is defined by peaks that share the same PeakContinuum pointer.
*/
struct RoiInfo
{
  std::vector<std::shared_ptr<const PeakDef>> peaks;  // sorted by mean
  double lower_energy;
  double upper_energy;
  size_t lower_channel;
  size_t upper_channel;
  double ref_energy;
  PeakContinuum::OffsetType offset_type;
  bool use_lls_for_cont;

  /** Whether the chi2 fit would solve this ROI's continuum by LLS; equal to `use_lls_for_cont`
   except for the all-in-Ceres objective, whose DOF and stamped Chi2DOF follow the chi2 fit's
   convention (see `dof_for_roi(...)`), so they mean the same whichever way the ROI is fit.
   */
  bool chi2_uses_lls_for_cont;
  double max_initial_sigma;   // max sigma across peaks in this ROI
  size_t num_fit_sigmas;      // number of peaks with fitFor(Sigma)==true

  /** Natural scale of a peak-CDF step coefficient for this ROI, in keV^-1.

   The Ceres variable is the dimensionless `step_coeff / cdf_step_scale`, which is the step's
   continuum-density change divided by the continuum level itself - O(0.01..1) for every ROI and
   detector, so a bound of +-2 means "the step may not move the continuum by more than twice its
   own level".  Without this the raw coefficient is O(1e-3..1e-6) and sits in the same parameter
   vector as means (~1) and continuum densities (~1e2), which invites premature termination
   against Ceres' parameter_tolerance.

   Computed once from the input peak amplitudes and the data, never from live fit parameters, so
   it is a constant and the Jet derivatives stay exact.  Left at 1.0 (a no-op) for every
   non-CDF-step ROI.
   */
  double cdf_step_scale = 1.0;

  /** Total counts in the ROI's channels; used to size a peak-amplitude parameter when the peak's
   own amplitude has not been fit yet (see `roi_amp_par_scale`).
   */
  double data_area = 0.0;
};//struct RoiInfo


/** Natural scale of the k-th peak-CDF step coefficient, in keV^-(1+k).

 `cdf_step_scale` is the scale of the constant term, which is in keV^-1.  Coefficient k multiplies
 `E'^k`, and `E'` spans the ROI, so its natural magnitude is smaller by `range^k` - without this,
 BiLinearStepCDF's keV^-2 slope would be scaled and bounded as though it were keV^-1, and the +-2
 bound would let it swing the step by `2*range` times the continuum level at the ROI's edge.
 */
double cdf_step_par_scale( const RoiInfo &roi, const size_t order )
{
  const double range = roi.upper_energy - roi.lower_energy;

  double scale = roi.cdf_step_scale;
  for( size_t i = 0; i < order; ++i )
    scale /= ((range > 0.0) ? range : 1.0);

  return scale;
}//double cdf_step_par_scale( const RoiInfo &, const size_t )


/** Natural scale of peak `i`'s amplitude parameter, so the Ceres variable is O(1). */
double roi_amp_par_scale( const RoiInfo &roi, const size_t peak_index )
{
  const double amp = roi.peaks[peak_index]->amplitude();
  if( std::isfinite(amp) && (amp > 1.0) )
    return amp;

  // A peak that has not been fit yet - or one that came back from the very defect this parameter
  //  block exists to fix, whose amplitude is exactly 1 - needs a scale from the data instead, or
  //  the parameter is in raw counts and its ceiling would sit absurdly low.
  const double sigma = roi.peaks[peak_index]->sigma();
  const double range = roi.upper_energy - roi.lower_energy;
  const double roi_counts = roi.data_area;
  const double frac = ((range > 0.0) && std::isfinite(sigma) && (sigma > 0.0))
                      ? (std::min)( 1.0, (2.0*sigma)/range ) : 1.0;

  const double est = frac * roi_counts;

  return (std::isfinite(est) && (est > 1.0)) ? est : 1.0;
}//double roi_amp_par_scale( const RoiInfo &, const size_t )


/** Number of peak-amplitude parameters this ROI contributes to the Ceres parameter vector.

 Normally zero: the amplitudes are linear in the model, so the least-squares solve inside
 `fit_amp_and_offset_imp(...)` recovers them for free and far more robustly than Ceres would.
 But that solve also solves the continuum polynomial, so it cannot be used when the caller has
 asked for a continuum coefficient to be held fixed - in that case the amplitudes have to be
 Ceres parameters like everything else, or they never get solved at all.
 */
size_t roi_amp_parameter_count( const RoiInfo &roi )
{
  if( roi.use_lls_for_cont )
    return 0;

  size_t nfit = 0;
  for( const std::shared_ptr<const PeakDef> &p : roi.peaks )
    nfit += p->fitFor( PeakDef::GaussAmplitude ) ? 1 : 0;

  return nfit;
}//size_t roi_amp_parameter_count( const RoiInfo &roi )


/** Number of continuum parameters this ROI contributes to the Ceres parameter vector.

 When the LLS solves the continuum, the polynomial terms are not Ceres parameters at all; only the
 peak-CDF step coefficients remain, because they are bilinear with the peak amplitudes and so
 cannot be part of the linear solve.  When it does not (any coefficient held fixed), every
 continuum parameter is a Ceres parameter.

 They occupy the first `roi_cont_parameter_count(roi)` slots of the ROI's parameter block.
 */
size_t roi_cont_parameter_count( const RoiInfo &roi )
{
  return roi.use_lls_for_cont ? PeakContinuum::num_cdf_step_pars( roi.offset_type )
                              : PeakContinuum::num_parameters( roi.offset_type );
}//size_t roi_cont_parameter_count( const RoiInfo &roi )


/** Computes `RoiInfo::cdf_step_scale` - the natural keV^-1 scale of a peak-CDF step coefficient.

 The step contributes `step_coeff * SUM_j(amp_j * CDF_j) * dx` counts to a channel, so
 `step_coeff * SUM(amp)` is a continuum count density.  Dividing the continuum's own density by
 the ROI peak area therefore gives a scale in the same units as `step_coeff`, and the ratio
 `step_coeff / scale` is the fractional change the step makes to the continuum.

 Both estimates deliberately fail toward being *too large*: overestimating the scale only loosens
 the parameter bound, whereas underestimating it would clamp a real step.  The result is always
 finite and strictly positive.
 */
double cdf_step_scale_for_roi( const RoiInfo &roi,
                               const std::shared_ptr<const SpecUtils::Measurement> &data )
{
  const double roi_width = roi.upper_energy - roi.lower_energy;
  if( !data || (roi.upper_channel < roi.lower_channel) || (roi_width <= 0.0) )
    return 1.0;

  const size_t nchan = 1 + roi.upper_channel - roi.lower_channel;

  // Continuum level from the ROI's edge channels - mostly continuum, little peak.
  const size_t nedge = std::max( size_t(1), std::min( size_t(4), nchan/4 ) );
  double low_sum = 0.0, up_sum = 0.0, low_width = 0.0, up_width = 0.0, roi_sum = 0.0;
  for( size_t i = 0; i < nchan; ++i )
  {
    const size_t channel = roi.lower_channel + i;
    const double counts = data->gamma_channel_content( channel );
    const double width = data->gamma_channel_width( channel );
    roi_sum += counts;

    if( i < nedge )
    {
      low_sum += counts;
      low_width += width;
    }

    if( (i + nedge) >= nchan )
    {
      up_sum += counts;
      up_width += width;
    }
  }//for( size_t i = 0; i < nchan; ++i )

  double density = 0.0;
  if( (low_width > 0.0) && (up_width > 0.0) )
    density = 0.5*((low_sum/low_width) + (up_sum/up_width));

  if( !std::isfinite(density) || (density <= 0.0) )
    density = roi_sum / roi_width;               // includes the peaks; an over-estimate

  if( !std::isfinite(density) || (density <= 0.0) )
    density = 1.0 / roi_width;                   // "one count across the ROI"; never zero

  // Total peak area in the ROI, including peaks whose amplitude is not being fit - the evaluator
  //  sums over all of them.
  double total_amp = 0.0;
  for( const std::shared_ptr<const PeakDef> &p : roi.peaks )
    total_amp += std::max( p->amplitude(), 0.0 );

  if( !std::isfinite(total_amp) || (total_amp < 1.0) )
  {
    // Fresh candidate peaks can still have zero amplitude; estimate it from the data so the
    //  resulting bound stays meaningful rather than becoming vacuous.
    total_amp = std::max( 1.0, roi_sum - (density * roi_width) );
  }

  const double scale = density / total_amp;

  return (!std::isfinite(scale) || (scale <= 0.0)) ? 1.0 : scale;
}//double cdf_step_scale_for_roi( const RoiInfo &, data )


/** PeakFitDiffCostFunction fits peaks from one or more ROIs simultaneously using Ceres.
 Equivalent of the MultiPeakFitChi2Fcn class, but using Levenberg-Marquardt differentiation.

 When multiple ROIs are present, skew parameters may be energy-dependent (interpolated between
 anchor energies at the spectrum lower/upper bounds), if the ROIs span more than 100 keV and
 the skew type has energy-dependent parameters.

 Parameter layout (default / shared-skew mode):
   [skew_lower_pars (num_skew values) | skew_upper_pars* (M values, only energy-dep params)]
   | ROI_0_cont | ROI_0_sigma | ROI_0_mean | ROI_0_amp** | ROI_1_cont | ... 
   * only present when m_fit_skew_energy_dependence == true

 Parameter layout (IndependentSkewValues option):
   ROI_0_cont | ROI_0_sigma | ROI_0_mean | ROI_0_amp** | ROI_0_skew (num_skew values)
   | ROI_1_cont | ROI_1_sigma | ROI_1_mean | ROI_1_amp** | ROI_1_skew | ...
   There is no shared skew block; each ROI carries its own num_skew parameters at the end of
   its per-ROI block.

   ** ROI_i_amp is present ONLY when that ROI's continuum is not being solved by the linear
      least-squares (i.e. a polynomial coefficient is pinned) - normally the amplitudes are
      recovered by the LLS and cost no Ceres parameters at all.  See roi_amp_parameter_count().
      Anything that indexes past the mean block must add it; see roi_skew_ptr / skew_base_idx.
*/
struct PeakFitDiffCostFunction
{
  /** A struct with the info to feed into Ceres. */
  struct ProblemSetup
  {
    /** The starting parameters to use. */
    vector<double> m_parameters;
    /** The indexes of paramters that need to be held constant in the fit.  E.g. if a mean is held constant for a peak whose FWHM/amp is being fit. */
    vector<int> m_constant_parameters;
    /** Lower bounds parameters should be allowed to go to; will not have a value if it shouldnt be restricted. Will be same size as `m_parameters`. */
    vector<std::optional<double>> m_lower_bounds;
    /** Upper bounds parameters should be allowed to go to; will not have a value if it shouldnt be restricted. Will be same size as `m_parameters`. */
    vector<std::optional<double>> m_upper_bounds;
  };//struct ProblemSetup


  /** Build the vector of RoiInfo structs from starting peaks and data.
   Groups peaks by shared PeakContinuum pointer, fills channel bounds, sigma, etc.
   Static so it can be called from the constructor initializer list via a lambda.
  */
  static std::vector<RoiInfo> make_rois(
      const std::shared_ptr<const SpecUtils::Measurement> &data,
      const std::vector<std::shared_ptr<const PeakDef>> &starting_peaks,
      const FitObjective objective )
  {
    // Group peaks by continuum pointer
    std::map<std::shared_ptr<const PeakContinuum>, std::vector<std::shared_ptr<const PeakDef>>> cont_to_peaks;
    for( const auto &p : starting_peaks )
      cont_to_peaks[p->continuum()].push_back( p );

    std::vector<RoiInfo> rois;
    rois.reserve( cont_to_peaks.size() );

    for( auto &kv : cont_to_peaks )
    {
      RoiInfo roi;
      roi.peaks = kv.second;
      std::sort( roi.peaks.begin(), roi.peaks.end(), &PeakDef::lessThanByMeanShrdPtr );

      const std::shared_ptr<const PeakContinuum> &cont = kv.first;
      roi.lower_energy = cont->lowerEnergy();
      roi.upper_energy = cont->upperEnergy();
      roi.lower_channel = data->find_gamma_channel( static_cast<float>( roi.lower_energy ) );
      roi.upper_channel = data->find_gamma_channel( static_cast<float>( roi.upper_energy ) );
      roi.offset_type = cont->type();

      const double prev_ref = cont->referenceEnergy();
      roi.ref_energy = ((prev_ref >= roi.lower_energy) && (prev_ref <= roi.upper_energy))
                       ? prev_ref : roi.lower_energy;

      // Use LLS for continuum if all continuum params are being fit (or NoOffset/External)
      if( (roi.offset_type == PeakContinuum::OffsetType::NoOffset)
          || (roi.offset_type == PeakContinuum::OffsetType::External) )
      {
        roi.use_lls_for_cont = true;
      }
      else
      {
        // Only a pinned *polynomial* coefficient forces us off the LLS path.  The peak-CDF step
        //  coefficients are Ceres parameters either way, so pinning one of those is honoured by
        //  simply marking it constant - see setup_roi_parameters(...).
        const size_t num_poly = PeakContinuum::num_linear_fit_pars( roi.offset_type );
        const vector<bool> cont_fit_for = cont->fitForParameter();

        roi.use_lls_for_cont = true;
        for( size_t i = 0; (i < num_poly) && (i < cont_fit_for.size()); ++i )
        {
          if( !cont_fit_for[i] )
          {
            roi.use_lls_for_cont = false;
            break;
          }
        }
      }

      // The all-in-Ceres objective makes every parameter a Ceres parameter (an External continuum
      //  is not one, so never gets that objective - see `fit_peaks_in_roi_imp(...)`).
      assert( (objective != FitObjective::PoissonAllInCeres)
              || (roi.offset_type != PeakContinuum::OffsetType::External) );
      roi.chi2_uses_lls_for_cont = roi.use_lls_for_cont;
      if( objective == FitObjective::PoissonAllInCeres )
        roi.use_lls_for_cont = false;


      // Max sigma across all peaks in this ROI
      roi.max_initial_sigma = 1.0;
      for( const auto &p : roi.peaks )
      {
        const double s = p->sigma();
        if( !isinf(s) && !isnan(s) && (s > 0.01) )
          roi.max_initial_sigma = std::max( roi.max_initial_sigma, s );
      }

      // Count how many peaks have Sigma being fit
      roi.num_fit_sigmas = 0;
      for( const auto &p : roi.peaks )
        roi.num_fit_sigmas += p->fitFor( PeakDef::Sigma ) ? 1 : 0;

      roi.data_area = 0.0;
      if( data && (roi.upper_channel >= roi.lower_channel) )
      {
        for( size_t ch = roi.lower_channel; ch <= roi.upper_channel; ++ch )
          roi.data_area += std::max( 0.0, static_cast<double>( data->gamma_channel_content(ch) ) );
      }

      // Scale that makes the peak-CDF step coefficient a dimensionless O(1) quantity; see the
      //  `RoiInfo::cdf_step_scale` documentation.
      if( PeakContinuum::num_cdf_step_pars( roi.offset_type ) )
        roi.cdf_step_scale = cdf_step_scale_for_roi( roi, data );

      rois.push_back( std::move(roi) );
    }//for( auto &kv : cont_to_peaks )

    // Sort ROIs by lower energy
    std::sort( rois.begin(), rois.end(), []( const RoiInfo &a, const RoiInfo &b ){
      return a.lower_energy < b.lower_energy;
    } );

    return rois;
  }//make_rois(...)


  PeakFitDiffCostFunction( const std::shared_ptr<const SpecUtils::Measurement> data,
                           const std::vector<std::shared_ptr<const PeakDef>> &starting_peaks,
                           const double roi_lower_energy,    // kept for API compat, ignored (ROI info comes from peaks)
                           const double roi_upper_energy,    // kept for API compat, ignored
                           const double continuum_ref_energy, // kept for API compat, ignored
                           const PeakDef::SkewType skew_type,
                           const PeakFitUtils::CoarseResolutionType det_type,
                           const Wt::WFlags<PeakFitLM::PeakFitLMOptions> options,
                           const FitObjective objective = FitObjective::NeymanChi2 )
  : m_data( data ),
    m_skew_type( skew_type ),
    m_det_type( det_type ),
    m_options( effective_options( options, data ) ),
    m_objective( objective ),
    m_ncalls( 0 ),
    m_rois( make_rois( data, starting_peaks, m_objective ) ),
    m_total_num_peaks( ([this]() -> size_t {
      size_t n = 0;
      for( const RoiInfo &r : m_rois )
        n += r.peaks.size();
      return n;
    })() ),
    m_fit_skew_energy_dependence( ([this, &skew_type]() -> bool {
      // IndependentSkewValues uses per-ROI skew blocks, not a shared energy-dependent block
      if( m_options.test( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
        return false;
      if( m_rois.size() <= 1 )
        return false;
      if( PeakDef::num_skew_parameters( skew_type ) == 0 )
        return false;
      // Check if any skew parameters are energy-dependent
      bool any_energy_dep = false;
      for( size_t i = 0; i < PeakDef::num_skew_parameters( skew_type ); ++i )
      {
        const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
        if( PeakDef::is_energy_dependent( skew_type, ct ) )
        {
          any_energy_dep = true;
          break;
        }
      }
      if( !any_energy_dep )
        return false;
      // Only energy-dependent if ROIs span more than 100 keV
      const double lower = m_rois.front().lower_energy;
      const double upper = m_rois.back().upper_energy;
      return (upper - lower) >= 100.0;
    })() ),
    m_skew_anchor_lower_energy( data->gamma_channel_lower( 0 ) ),
    m_skew_anchor_upper_energy( data->gamma_channel_upper( data->num_gamma_channels() - 1 ) ),
    m_external_continuum( ([this]() -> std::shared_ptr<const SpecUtils::Measurement> {
      for( const RoiInfo &roi : m_rois )
      {
        if( roi.offset_type == PeakContinuum::OffsetType::External )
        {
          for( const auto &p : roi.peaks )
          {
            if( p->continuum()->externalContinuum() )
              return p->continuum()->externalContinuum();
          }
        }
      }
      return nullptr;
    })() ),
    m_num_parameters( number_parameters() ),
    m_num_residuals( number_residuals() )
#if( PEAK_FIT_LM_PARALLEL_ROIS )
    , m_thread_pool( (m_rois.size() > 2)
        ? std::make_unique<boost::asio::thread_pool>(
            std::min( m_rois.size(),
                      static_cast<size_t>( std::max( 1, SpecUtilsAsync::num_physical_cpu_cores() ) ) ) )
        : nullptr )
#endif
  {
    if( !m_data || !m_data->gamma_counts() || m_data->gamma_counts()->empty() )
      throw runtime_error( "PeakFitDiffCostFunction: !m_data or !gamma_counts()" );

    if( m_rois.empty() )
      throw runtime_error( "PeakFitDiffCostFunction: no ROIs to fit" );

    const std::vector<float> &counts = *m_data->gamma_counts();
    const std::shared_ptr<const vector<float>> &energies = m_data->channel_energies();

    for( const RoiInfo &roi : m_rois )
    {
      if( roi.lower_energy >= roi.upper_energy )
        throw runtime_error( "PeakFitDiffCostFunction: ROI lower_energy >= upper_energy" );
      if( roi.upper_channel < roi.lower_channel )
        throw runtime_error( "PeakFitDiffCostFunction: ROI upper_channel < lower_channel" );
      if( roi.upper_channel >= counts.size() )
        throw runtime_error( "PeakFitDiffCostFunction: ROI upper_channel >= counts.size()" );
      if( !energies || (roi.upper_channel >= energies->size()) )
        throw runtime_error( "PeakFitDiffCostFunction: ROI upper_channel >= energies.size()" );
      if( roi.peaks.empty() )
        throw runtime_error( "PeakFitDiffCostFunction: ROI has no peaks" );
      if( roi.peaks.size() >= (roi.upper_channel - roi.lower_channel) )
        throw runtime_error( "PeakFitDiffCostFunction: too many peaks for ROI channel range" );

      if( (roi.offset_type == PeakContinuum::OffsetType::External)
          && (!m_external_continuum
              || (m_external_continuum->num_gamma_channels() < 7)
              || !m_external_continuum->energy_calibration()
              || !m_external_continuum->energy_calibration()->valid()) )
        throw runtime_error( "PeakFitDiffCostFunction: external continuum wanted but not valid" );
    }//for( const RoiInfo &roi : m_rois )

    // Check all ROIs have positive DOF
    for( size_t i = 0; i < m_rois.size(); ++i )
    {
      if( dof_for_roi( i ) <= 0.0 )
        throw runtime_error( "PeakFitDiffCostFunction: not enough data channels to fit parameters in ROI" );
    }
  }//PeakFitDiffCostFunction constructor


  // Returns count of sigma parameters for one ROI
  size_t roi_sigma_parameter_count( const RoiInfo &roi ) const
  {
    if( m_options.test( PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent ) )
      return roi.num_fit_sigmas;
    else
      return std::min( roi.num_fit_sigmas, size_t(2) );
  }

  // Returns total parameter count for one ROI (cont + sigma + mean params, plus per-ROI skew
  // when IndependentSkewValues is set).
  size_t roi_parameter_count( const RoiInfo &roi ) const
  {
    const size_t cont_pars = roi_cont_parameter_count( roi );
    size_t n = cont_pars + roi_sigma_parameter_count( roi ) + roi.peaks.size()
               + roi_amp_parameter_count( roi );
    if( m_options.test( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
      n += PeakDef::num_skew_parameters( m_skew_type );
    return n;
  }

  // Returns total residual count for one ROI (data channels + punishment residuals)
  size_t roi_residual_count( const RoiInfo &roi ) const
  {
    const size_t n_chan = roi.upper_channel - roi.lower_channel + 1;
    const bool punish_close = !m_options.test( PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose );
#if( ENABLE_PUNISH_STAT_INSIG_PEAKS )
    const bool punish_insig = m_options.test( PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig );
#else
    const bool punish_insig = false;
#endif
    size_t n = n_chan;
    if( (punish_close || punish_insig) && roi.peaks.size() > 1 )
      n += roi.peaks.size() - 1;
    if( punish_insig )
      n += 1;
    return n;
  }

  // Returns count of shared skew parameters (doubles for energy-dep params when multi-ROI).
  // Returns 0 when IndependentSkewValues is set, since skew params are folded into each ROI block.
  size_t skew_parameter_count() const
  {
    if( m_options.test( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
      return 0;

    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );
    if( !m_fit_skew_energy_dependence )
      return num_skew;

    // Layout: [sp0..spN-1 (lower/base values)] [upper values for energy-dep params only]
    size_t num_energy_dep = 0;
    for( size_t i = 0; i < num_skew; ++i )
    {
      const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
      if( PeakDef::is_energy_dependent( m_skew_type, ct ) )
        num_energy_dep += 1;
    }
    return num_skew + num_energy_dep;
  }

  // Returns degrees of freedom for a single ROI.
  // Shared skew parameters are not counted here (they span all ROIs).
  // When IndependentSkewValues is set, per-ROI skew parameters are counted here.
  //
  // KNOWN DISCREPANCY, pre-existing and deliberately left alone (2026-09): on the LLS path
  //  `roi_cont_parameter_count(roi)` returns only `num_cdf_step_pars(type)`, which is zero for
  //  every non-CDF continuum type, so the 1-4 polynomial coefficients that
  //  `PeakFit::fit_amp_and_offset_imp(...)` solves never reduce DOF.  They are fit, so they do
  //  consume it, and this overstates DOF by `num_linear_fit_pars(type)` for every ROI on this
  //  path (it is decided by `chi2_uses_lls_for_cont`, so the all-in-Ceres objective stamps the
  //  same convention).  `get_chi2_and_dof_for_roi(...)` in src/PeakFit.cpp counts every free continuum
  //  parameter and is the correct convention; see the long note there for why the fix has not been
  //  applied (the stamped chi2/DOF gates automated peak acceptance) and for the other, deliberate,
  //  divergences between the two.  If you fix this, count `continuum()->fitForParameter()`
  //  directly - do NOT change `roi_cont_parameter_count()`, which also sizes and indexes the Ceres
  //  parameter blocks.
  double dof_for_roi( const size_t roi_index ) const
  {
    assert( roi_index < m_rois.size() );
    const RoiInfo &roi = m_rois[roi_index];
    const double num_channels = static_cast<double>( 1 + roi.upper_channel - roi.lower_channel );
    const size_t num_fit_cont = roi.chi2_uses_lls_for_cont ? PeakContinuum::num_cdf_step_pars( roi.offset_type )
                                                           : PeakContinuum::num_parameters( roi.offset_type );

    size_t num_fixed = 0;
    for( const auto &p : roi.peaks )
    {
      num_fixed += p->fitFor( PeakDef::GaussAmplitude ) ? 0 : 1;
      num_fixed += p->fitFor( PeakDef::Mean ) ? 0 : 1;
    }
    {
      // Continuum coefficients the caller pinned do not consume a degree of freedom.  On the LLS
      //  path only the peak-CDF step coefficients are ours to count (the polynomial terms are not
      //  in `num_fit_cont` there); off it, every continuum parameter is.
      const vector<bool> cont_fit_for = roi.peaks[0]->continuum()->fitForParameter();
      const size_t first = roi.chi2_uses_lls_for_cont
                           ? PeakContinuum::num_linear_fit_pars( roi.offset_type ) : size_t(0);
      for( size_t i = first; i < cont_fit_for.size(); ++i )
        num_fixed += cont_fit_for[i] ? 0 : 1;
    }

    // Count fitted per-ROI skew parameters when IndependentSkewValues is active.
    // Mirrors setup_roi_parameters exactly: a parameter is counted as fitted if ANY
    // matching-type peak in the ROI has fitFor==true for it (OR semantics).
    size_t num_fit_skew = 0;
    if( m_options.test( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
    {
      const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skew_type );
      if( num_skew_pars > 0 )
      {
        vector<bool> fit_skew( num_skew_pars, false );
        for( const auto &p : roi.peaks )
        {
          if( p->skewType() != m_skew_type )
            continue;
          for( size_t i = 0; i < num_skew_pars; ++i )
          {
            const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
            if( p->fitFor( ct ) )
              fit_skew[i] = true;
          }
        }
        for( const bool fit : fit_skew )
          num_fit_skew += fit ? 1 : 0;
      }
    }

    return num_channels
           - 2.0*static_cast<double>( roi.peaks.size() )
           + static_cast<double>( num_fixed )
           - static_cast<double>( roi_sigma_parameter_count( roi ) )
           - static_cast<double>( num_fit_cont )
           - static_cast<double>( num_fit_skew );
  }

  // Returns all starting peaks flattened across all ROIs (for skew param initialization)
  std::vector<std::shared_ptr<const PeakDef>> all_starting_peaks() const
  {
    std::vector<std::shared_ptr<const PeakDef>> all;
    for( const RoiInfo &roi : m_rois )
      for( const auto &p : roi.peaks )
        all.push_back( p );
    return all;
  }

  // Computes total parameter count from scratch (called from initializer list)
  size_t number_parameters() const
  {
    size_t n = skew_parameter_count();
    for( const RoiInfo &roi : m_rois )
      n += roi_parameter_count( roi );
    return n;
  }

  // Computes total residual count from scratch (called from initializer list)
  size_t number_residuals() const
  {
    size_t n = 0;
    for( const RoiInfo &roi : m_rois )
      n += roi_residual_count( roi );
    return n;
  }


  //
  template<typename T>
  static std::shared_ptr<T> create_continuum( const std::shared_ptr<T> &other_cont [[maybe_unused]]) {
      return std::make_shared<T>();
  }

  /** Apply the shared skew parameters to all peaks in the given vector.
   For single-ROI (or no energy dependence), all peaks get the same skew values.
   For multi-ROI with energy-dependent params, the skew is interpolated between
   the lower and upper anchor energies based on each peak's mean energy.

   Skew parameter layout in params[0..skew_block_size-1]:
     params[0..num_skew-1]:        base (lower-anchor) values for all skew params
     params[num_skew..num_skew+M-1]: upper-anchor values for energy-dep params only (M = count of energy-dep params)

   uncertainties: if non-null, has the same layout as params (diagonal of covariance); used to
     set skew uncertainty on each peak for non-energy-dependent params.
   covariance: if non-null, the full num_total_pars x num_total_pars row-major covariance matrix
     for the entire parameter vector; needed for proper error propagation of interpolated
     energy-dependent skew params.  The indices in this matrix correspond to the global parameter
     array, so params[0] corresponds to covariance[0][0], etc.
   params_offset: the index of params[0] in the global parameter array (needed to index covariance).
  */
  template<typename PeakType, typename T>
  void apply_skew_to_peaks( vector<PeakType> &peaks, const T * const params,
                            const T * const uncertainties = nullptr,
                            const double * const covariance = nullptr,
                            const size_t num_total_pars = 0,
                            const size_t params_offset = 0 ) const
  {
    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );
    if( num_skew == 0 )
      return;

    const double energy_span = m_skew_anchor_upper_energy - m_skew_anchor_lower_energy;

    for( PeakType &peak : peaks )
    {
      peak.setSkewType( m_skew_type );

      if( !m_fit_skew_energy_dependence )
      {
        // Single ROI or no energy dependence: all peaks share the same skew values
        for( size_t i = 0; i < num_skew; ++i )
        {
          const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
          peak.set_coefficient( params[i], ct );
          if( uncertainties )
            peak.set_uncertainty( uncertainties[i], ct );
        }
      }
      else
      {
        // Multi-ROI with energy dependence: interpolate per-peak based on mean
        T mean_val;
        if constexpr ( std::is_same_v<T, double> )
          mean_val = T( peak.mean() );
        else
          mean_val = peak.mean();

        const T mean_frac = (mean_val - T(m_skew_anchor_lower_energy)) / T(energy_span);

        size_t upper_offset = num_skew; // index of first upper-anchor param
        for( size_t i = 0; i < num_skew; ++i )
        {
          const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
          T val;
          if( PeakDef::is_energy_dependent( m_skew_type, ct ) )
          {
            // Interpolate between lower (params[i]) and upper (params[upper_offset])
            val = params[i] + mean_frac * (params[upper_offset] - params[i]);
            peak.set_coefficient( val, ct );

            // Uncertainty propagation only applies when T=double (post-solve), not during Jet-based solve.
            if constexpr ( std::is_same_v<T, double> )
            {
              // Proper error propagation: val = (1-f)*p_lower + f*p_upper
              // sigma^2 = (1-f)^2*Var[lower] + f^2*Var[upper] + 2*(1-f)*f*Cov[lower,upper]
              if( covariance && (num_total_pars > 0) )
              {
                const double f = mean_frac;
                const double one_minus_f = 1.0 - f;
                const size_t gi = params_offset + i;              // global index of lower-anchor param
                const size_t gu = params_offset + upper_offset;   // global index of upper-anchor param
                const double var_lower = covariance[gi * num_total_pars + gi];
                const double var_upper = covariance[gu * num_total_pars + gu];
                const double cov_lu    = covariance[gi * num_total_pars + gu];
                const double variance  = one_minus_f*one_minus_f*var_lower
                                         + f*f*var_upper
                                         + 2.0*one_minus_f*f*cov_lu;
                if( variance > 0.0 )
                  peak.set_uncertainty( sqrt(variance), ct );
              }
              else if( uncertainties )
              {
                // Fallback: propagate in quadrature ignoring correlation
                const double f = mean_frac;
                const double sigma_lower = uncertainties[i];
                const double sigma_upper = uncertainties[upper_offset];
                const double sigma = sqrt( (1.0-f)*(1.0-f)*sigma_lower*sigma_lower
                                           + f*f*sigma_upper*sigma_upper );
                if( sigma > 0.0 )
                  peak.set_uncertainty( sigma, ct );
              }
            }//if constexpr T==double

            upper_offset += 1;
          }
          else
          {
            val = params[i];
            peak.set_coefficient( val, ct );
            if( uncertainties )
              peak.set_uncertainty( uncertainties[i], ct );
          }
        }
      }
    }//for( PeakType &peak : peaks )
  }//apply_skew_to_peaks(...)


  template<typename PeakType,typename T>
  vector<PeakType> parametersToPeaks( const T * const params, const T * const uncertainties, T *residuals,
                                      const double * const covariance = nullptr,
                                      const size_t num_total_pars = 0,
                                      std::vector<std::vector<double>> * const roi_models = nullptr ) const
  {
    // Each ROI writes only its own entry, so the ROIs may be processed concurrently.
    if( roi_models )
      roi_models->assign( m_rois.size(), std::vector<double>() );

    std::unique_ptr<vector<T>> local_residuals;
    if( !residuals )
    {
      local_residuals.reset( new vector<T>( number_residuals(), T(0.0) ) );
      residuals = local_residuals->data();
    }
    else
    {
      const size_t nresids = number_residuals();
      for( size_t i = 0; i < nresids; ++i )
        residuals[i] = T(0.0);
    }

    const size_t num_skew = PeakDef::num_skew_parameters( m_skew_type );
    const size_t skew_block_size = skew_parameter_count();

    // Pre-compute per-ROI parameter and residual offsets so we can launch ROIs in parallel.
    const size_t nrois = m_rois.size();
    vector<size_t> roi_param_offsets( nrois ), roi_residual_offsets( nrois );
    {
      size_t poff = skew_block_size, roff = 0;
      for( size_t i = 0; i < nrois; ++i )
      {
        roi_param_offsets[i]    = poff;
        roi_residual_offsets[i] = roff;
        poff += roi_parameter_count( m_rois[i] );
        roff += roi_residual_count( m_rois[i] );
      }
    }

    // Lambda that processes one ROI and returns its peaks; writes residuals into roi_residuals.
    // All inputs are read-only except roi_residuals (which is per-ROI private storage).
    // params[] is read-only (shared skew block + ROI params), safe to read from multiple threads.
    const auto process_one_roi = [&]( const size_t roi_idx,
                                      const size_t param_offset,
                                      T * const roi_residuals ) -> vector<PeakType>
    {
      const RoiInfo &roi = m_rois[roi_idx];
      const size_t num_roi_peaks = roi.peaks.size();
      const size_t num_sigmas_fit = roi_sigma_parameter_count( roi );
      const size_t num_fit_cont = roi_cont_parameter_count( roi );
      const size_t num_amps_fit = roi_amp_parameter_count( roi );
      const size_t nchannel = roi.upper_channel - roi.lower_channel + 1;
      const double range = roi.upper_energy - roi.lower_energy;

      // Pointer into the parameter array for this ROI
      const T * const roi_params = params + param_offset;

      vector<PeakType> peaks, fixed_amp_peaks;

      // Compute min/max means across this ROI for sigma interpolation
      T min_mean = T( std::numeric_limits<double>::infinity() );
      T max_mean = T( -std::numeric_limits<double>::infinity() );
      for( size_t i = 0; i < num_roi_peaks; ++i )
      {
        const size_t mean_par_index = num_fit_cont + num_sigmas_fit + i;
        const T frac = roi_params[mean_par_index] - T(0.5);
        const T mean = T(roi.lower_energy) + frac * T(range);
        min_mean = min( min_mean, mean );
        max_mean = max( max_mean, mean );
      }

      size_t fit_sigma_num = 0, fit_amp_num = 0;
      for( size_t i = 0; i < num_roi_peaks; ++i )
      {
        const std::shared_ptr<const PeakDef> &src_peak = roi.peaks[i];
        const size_t mean_par_index = num_fit_cont + num_sigmas_fit + i;
        assert( mean_par_index < roi_parameter_count( roi ) );

        const bool fit_amp = src_peak->fitFor( PeakDef::GaussAmplitude );

        const T frac = roi_params[mean_par_index] - T(0.5);
        const T mean = T(roi.lower_energy) + frac * T(range);

        // With the LLS in play the amplitude is a placeholder it will solve for; without it, the
        //  amplitude is a Ceres parameter of its own (scaled so the variable is O(1)).
        T amp = T(1.0), amp_uncert = T(0.0);
        if( !fit_amp )
        {
          amp = T( src_peak->amplitude() );
        }else if( num_amps_fit )
        {
          const size_t amp_index = num_fit_cont + num_sigmas_fit + num_roi_peaks + fit_amp_num;
          assert( amp_index < roi_parameter_count( roi ) );

          const double amp_scale = roi_amp_par_scale( roi, i );
          amp = roi_params[amp_index] * T(amp_scale);
          if( uncertainties )
            amp_uncert = uncertainties[param_offset + amp_index] * T(amp_scale);
        }

        T sigma, sigma_uncert;

        if( !src_peak->fitFor( PeakDef::Sigma ) || (num_sigmas_fit == 0) )
        {
          sigma = T( src_peak->sigma() );
          sigma_uncert = T( src_peak->sigmaUncert() );
        }
        else if( m_options.test( PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent ) )
        {
          assert( fit_sigma_num < roi.num_fit_sigmas );
          const size_t sigma_index = num_fit_cont + fit_sigma_num;
          assert( sigma_index < roi_parameter_count( roi ) );

          sigma = roi_params[sigma_index] * T(roi.max_initial_sigma);
          sigma_uncert = uncertainties
                         ? (uncertainties[param_offset + sigma_index] * T(roi.max_initial_sigma))
                         : T(0.0);

          if( isnan(sigma) || isinf(sigma) )
            throw runtime_error( "Inf or NaN sigma (AllPeakFwhmIndependent)" );
        }
        else
        {
          // First sigma param = absolute sigma for lowest-energy peak;
          // second sigma param (if present) = relative multiplier for highest-energy peak.
          const size_t sigma_index = num_fit_cont;
          sigma = roi_params[sigma_index] * T(roi.max_initial_sigma);
          sigma_uncert = uncertainties
                         ? (uncertainties[param_offset + sigma_index] * T(roi.max_initial_sigma))
                         : T(0.0);

          if( (i > 0) && (num_sigmas_fit > 1) )
          {
            T frac_dist;
            const T mean_dists = max_mean - min_mean;
            if( mean_dists > T(1.0E-6) )
              frac_dist = (mean - min_mean) / mean_dists;
            else
              frac_dist = T(0.5);

            const T first_sigma = sigma;
            const T last_sigma  = sigma * roi_params[sigma_index + 1];
            sigma = first_sigma + frac_dist * (last_sigma - first_sigma);

            if( uncertainties )
            {
              const T first_uncert = sigma_uncert;
              const T last_uncert  = sigma_uncert * uncertainties[param_offset + sigma_index + 1];
              sigma_uncert = first_uncert + frac_dist * (last_uncert - first_uncert);
            }
          }//if( (i > 0) && (num_sigmas_fit > 1) )

          if( isnan(sigma) || isinf(sigma) )
            throw runtime_error( "Inf or NaN sigma (shared sigma)" );
        }//sigma selection

        // Protect against sigma smaller than one channel width (only when fitting sigma)
        if( (sigma < 0.2) && src_peak->fitFor( PeakDef::Sigma ) && (num_sigmas_fit != 0) )
        {
          const float left_val  = m_data->gamma_channel_lower( roi.lower_channel );
          const float upper_val = m_data->gamma_channel_upper( roi.upper_channel );
          const double avrg_chnl_width = (upper_val - left_val) / static_cast<double>( nchannel );
          const double min_reasonable = (0.42*0.99) * avrg_chnl_width;

          double sigma_d;
          if constexpr ( std::is_same_v<T, double> )
            sigma_d = sigma;
          else
            sigma_d = sigma.a;

          if( (sigma_d <= min_reasonable) || (sigma_d < 0.00025) )
            throw runtime_error( "peak sigma (" + std::to_string(sigma_d)
                                 + ") smaller than reasonable ("
                                 + std::to_string( min_reasonable ) + ")" );
        }

        if( src_peak->fitFor( PeakDef::Sigma ) )
          fit_sigma_num += 1;
        if( fit_amp )
          fit_amp_num += 1;

        PeakType peak;
        peak.setMean( mean );
        peak.setSigma( sigma );
        peak.setAmplitude( amp );
        // Skew applied later via apply_skew_to_peaks
        peak.setSkewType( m_skew_type );

        if( uncertainties )
        {
          T mean_uncert = uncertainties[param_offset + mean_par_index];
          if( !src_peak->fitFor( PeakDef::Mean ) )
            mean_uncert = T( src_peak->meanUncert() );
          if( mean_uncert > 0.0 )
            peak.setMeanUncert( mean_uncert );
          if( sigma_uncert > 0.0 )
            peak.setSigmaUncert( sigma_uncert );
        }

        if( fit_amp )
        {
          if( num_amps_fit && (amp_uncert > 0.0) )
            peak.setAmplitudeUncert( amp_uncert );
          peaks.push_back( peak );
        }
        else
        {
          peak.setAmplitudeUncert( T( src_peak->amplitudeUncert() ) );
          fixed_amp_peaks.push_back( peak );
        }
      }//for( size_t i = 0; i < num_roi_peaks; ++i )

      assert( fit_sigma_num == roi.num_fit_sigmas );
      std::sort( begin(peaks), end(peaks), &PeakType::lessThanByMean );

      // In IndependentSkewValues mode the skew params are at the tail of this ROI's own block;
      // otherwise they live at params[0..num_skew-1] (the shared / energy-dependent block).
      const T *roi_skew_ptr;
      if( m_options.test( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
        roi_skew_ptr = roi_params + num_fit_cont + num_sigmas_fit + num_roi_peaks + num_amps_fit;
      else
        roi_skew_ptr = params;

      // The skew must be on the peaks before any model is computed from them: the non-LLS branch
      //  below, and the fixed-amplitude peaks in the LLS, evaluate `gauss_integral(...)` directly
      //  (without it a skewed peak is evaluated with unset skew parameters - NaN).  Applied again,
      //  with uncertainties, once the peaks are complete.
      apply_skew_to_peaks<PeakType, T>( peaks, roi_skew_ptr );
      apply_skew_to_peaks<PeakType, T>( fixed_amp_peaks, roi_skew_ptr );

      // --- Compute predicted channel counts for this ROI ---
      const shared_ptr<const vector<float>> &energies_ptr = m_data->channel_energies();
      const vector<float> &counts_vec = *m_data->gamma_counts();
      const float * const energies       = energies_ptr->data() + roi.lower_channel;
      const float * const channel_counts = counts_vec.data()    + roi.lower_channel;

      assert( !peaks.empty() || !fixed_amp_peaks.empty() );
      auto continuum = create_continuum( (!peaks.empty())
                                         ? peaks.front().getContinuum()
                                         : fixed_amp_peaks.front().getContinuum() );
      continuum->setRange( T(roi.lower_energy), T(roi.upper_energy) );
      continuum->setType( roi.offset_type );

      vector<T> peak_counts( nchannel, T(0.0) );

      // The skew parameter vector to pass to fit_amp_and_offset_imp.
      const vector<T> skew_pars( roi_skew_ptr, roi_skew_ptr + num_skew );

      // For a ROI being refit by IRLS, each channel's variance (see `m_irls_variances`); null for the
      //  chi2 fit, including the first pass of the IRLS fit.
      const bool irls_roi = (roi_idx < m_irls_variances.size()) && !m_irls_variances[roi_idx].empty();
      assert( !irls_roi || (m_irls_variances[roi_idx].size() == nchannel) );
      const double * const irls_variances = irls_roi ? m_irls_variances[roi_idx].data() : nullptr;

      if( roi.offset_type == PeakContinuum::OffsetType::External )
      {
        assert( m_external_continuum );
        continuum->setExternalContinuum( m_external_continuum );

        if( !peaks.empty() && irls_variances )
        {
          // The external continuum is a fixed part of the model; nothing is subtracted or clipped.
          vector<double> external_counts( nchannel );
          for( size_t i = 0; i < nchannel; ++i )
            external_counts[i] = m_external_continuum->gamma_integral( energies[i], energies[i+1] );

          vector<T> means, sigmas;
          for( const auto &p : peaks )
          {
            means.push_back( p.mean() );
            sigmas.push_back( p.sigma() );
          }

          vector<T> amplitudes, cont_coeffs, amp_uncerts, cont_uncerts;
          PeakFitLMObjective::fit_amp_and_offset_weighted<PeakType,T>( energies, channel_counts,
                                      irls_variances, external_counts.data(), nchannel,
                                      roi.offset_type, nullptr, T(roi.ref_energy),
                                      means, sigmas, fixed_amp_peaks, m_skew_type, skew_pars.data(),
                                      amplitudes, cont_coeffs, amp_uncerts, cont_uncerts, &peak_counts[0] );

          assert( peaks.size() == amplitudes.size() );
          for( size_t pi = 0; pi < peaks.size(); ++pi )
          {
            peaks[pi].setAmplitude( amplitudes[pi] );
            peaks[pi].setAmplitudeUncert( amp_uncerts[pi] );
          }
        }
        else if( !peaks.empty() )
        {
          // Subtract external continuum from data for amplitude fitting
          vector<float> data_copy( channel_counts, channel_counts + nchannel );
          vector<float> data_variances( channel_counts, channel_counts + nchannel );
          for( size_t i = 0; i < nchannel; ++i )
          {
            data_copy[i] -= m_external_continuum->gamma_integral( energies[i], energies[i+1] );
            data_copy[i] = std::max( 0.0f, data_copy[i] );
            if( data_variances[i] < PEAK_FIT_MIN_CHANNEL_UNCERT )
              data_variances[i] = PEAK_FIT_MIN_CHANNEL_UNCERT;
          }

          vector<T> means, sigmas;
          for( const auto &p : peaks )
          {
            means.push_back( p.mean() );
            sigmas.push_back( p.sigma() );
          }

          vector<T> amplitudes, cont_coeffs, amp_uncerts, cont_uncerts;
          {
            // An External continuum has no CDF step coefficients, so there is nothing to pass.
            assert( !PeakContinuum::num_cdf_step_pars( roi.offset_type ) );
            PeakFit::fit_amp_and_offset_imp( energies, &data_copy[0], &data_variances[0], nchannel,
                                             roi.offset_type, nullptr, T(roi.ref_energy),
                                             means, sigmas, fixed_amp_peaks, m_skew_type, skew_pars.data(),
                                             amplitudes, cont_coeffs, amp_uncerts, cont_uncerts,
                                             &peak_counts[0] );
          }

          assert( peaks.size() == amplitudes.size() );
          for( size_t pi = 0; pi < peaks.size(); ++pi )
          {
            peaks[pi].setAmplitude( amplitudes[pi] );
            if( !amp_uncerts.empty() )
              peaks[pi].setAmplitudeUncert( amp_uncerts[pi] );
          }
        }
        else
        {
          for( PeakType &fp : fixed_amp_peaks )
            fp.gauss_integral( energies, &peak_counts[0], nchannel );
        }

        if( peaks.empty() || !irls_variances )  //else the model already includes it
        {
          for( size_t i = 0; i < nchannel; ++i )
            peak_counts[i] += T( m_external_continuum->gamma_integral( energies[i], energies[i+1] ) );
        }
      }
      else if( !roi.use_lls_for_cont )
      {
        // Non-LLS continuum: continuum parameters are in roi_params[0..num_fit_cont-1].  Their
        //  uncertainties are only passed on for the all-in-Ceres objective, leaving the default
        //  fit's behaviour unchanged.
        //  (A NoOffset continuum - on this branch only for that objective - has no parameters.)
        const bool set_cont_uncerts = uncertainties && (m_objective == FitObjective::PoissonAllInCeres);
        if( num_fit_cont )
          continuum->setParameters( T(roi.ref_energy), roi_params,
                                    set_cont_uncerts ? (uncertainties + param_offset) : nullptr );

        for( PeakType &p : peaks )
          p.gauss_integral( energies, &peak_counts[0], nchannel );
        for( PeakType &fp : fixed_amp_peaks )
          fp.gauss_integral( energies, &peak_counts[0], nchannel );

        if( PeakContinuum::is_peak_cdf_step_continuum( roi.offset_type ) )
        {
          // A CDF step continuum's step term is defined by the ROI's peaks, so it needs the
          //  peak-aware integrator; the peaks-less one throws for these types.
          vector<const PeakType *> all_peaks;
          all_peaks.reserve( peaks.size() + fixed_amp_peaks.size() );
          for( const PeakType &p : peaks )
            all_peaks.push_back( &p );
          for( const PeakType &fp : fixed_amp_peaks )
            all_peaks.push_back( &fp );

          PeakDists::offset_integral( *continuum, energies, &peak_counts[0], nchannel, m_data,
                                      all_peaks.data(), all_peaks.size() );
        }else
        {
          PeakDists::offset_integral( *continuum, energies, &peak_counts[0], nchannel, m_data );
        }
      }
      else
      {
        // LLS continuum: fit continuum coefficients simultaneously with amplitudes
        vector<T> means, sigmas;
        for( const auto &p : peaks )
        {
          means.push_back( p.mean() );
          sigmas.push_back( p.sigma() );
        }

        // Ceres holds the step coefficients in dimensionless form; convert each back to its own
        //  physical units (keV^-1 for the constant term, keV^-2 for the energy slope).
        vector<T> step_coeffs( num_fit_cont );
        for( size_t k = 0; k < num_fit_cont; ++k )
          step_coeffs[k] = roi_params[k] * T(cdf_step_par_scale( roi, k ));

        vector<T> amplitudes, cont_coeffs, amp_uncerts, cont_uncerts;
        if( irls_variances )
        {
          PeakFitLMObjective::fit_amp_and_offset_weighted<PeakType,T>( energies, channel_counts,
                                        irls_variances, nullptr, nchannel,
                                        roi.offset_type, step_coeffs.empty() ? nullptr : step_coeffs.data(),
                                        T(roi.ref_energy), means, sigmas, fixed_amp_peaks, m_skew_type,
                                        skew_pars.data(), amplitudes, cont_coeffs, amp_uncerts, cont_uncerts,
                                        &peak_counts[0] );
        }else
        {
          PeakFit::fit_amp_and_offset_imp( energies, channel_counts, nullptr, nchannel,
                                           roi.offset_type,
                                           step_coeffs.empty() ? nullptr : step_coeffs.data(),
                                           T(roi.ref_energy),
                                           means, sigmas, fixed_amp_peaks, m_skew_type, skew_pars.data(),
                                           amplitudes, cont_coeffs, amp_uncerts, cont_uncerts,
                                           &peak_counts[0] );
        }
        assert( peaks.size() == amplitudes.size() );
        for( size_t pi = 0; pi < peaks.size(); ++pi )
        {
          peaks[pi].setAmplitude( amplitudes[pi] );
          if( !amp_uncerts.empty() )
            peaks[pi].setAmplitudeUncert( amp_uncerts[pi] );
        }

        if( (roi.offset_type != PeakContinuum::OffsetType::NoOffset)
            && (roi.offset_type != PeakContinuum::OffsetType::External) )
        {
          // fit_amp_and_offset_imp returns only the polynomial coefficients, but the continuum
          //  expects the step coefficients appended.  The rescale is linear, so the uncertainties
          //  scale by exactly the same factor.
          for( size_t k = 0; k < num_fit_cont; ++k )
          {
            cont_coeffs.push_back( step_coeffs[k] );
            cont_uncerts.push_back( uncertainties
                        ? (T(uncertainties[param_offset + k]) * T(cdf_step_par_scale( roi, k )))
                        : T(0.0) );
          }

          continuum->setParameters( T(roi.ref_energy), cont_coeffs.data(), cont_uncerts.data() );
        }
      }//continuum handling

      // Merge fixed-amp peaks back into peaks
      if( !fixed_amp_peaks.empty() )
      {
        for( PeakType &fp : fixed_amp_peaks )
          peaks.push_back( fp );
        std::sort( begin(peaks), end(peaks), &PeakType::lessThanByMean );
      }

      // Set continuum and apply skew parameters to all peaks in this ROI.
      // In IndependentSkewValues mode, pass the per-ROI skew pointer; otherwise pass params
      // (which points to the shared/energy-dependent skew block at the start of the global array).
      for( PeakType &p : peaks )
        p.setContinuum( continuum );

      // Pass uncertainties at the same offset as the skew params, so uncertainties are set too.
      const T *roi_skew_uncert_ptr = uncertainties ? (uncertainties + (roi_skew_ptr - params)) : nullptr;
      const size_t skew_params_offset = static_cast<size_t>( roi_skew_ptr - params );
      apply_skew_to_peaks<PeakType, T>( peaks, roi_skew_ptr, roi_skew_uncert_ptr,
                                        covariance, num_total_pars, skew_params_offset );

      // --- Compute residuals for this ROI ---
      const double deviance_floor = sm_deviance_min_expected_frac
                                   * std::max( 1.0, roi.data_area / static_cast<double>(nchannel) );

      // `chi2` is always the modified-Neyman chi2, which is what gets stamped as Chi2DOF, however
      //  the ROI is fit; see the notes on the fit statistic in `PeakFitLMOptions`.
      T chi2( 0.0 );
      for( size_t ch = 0; ch < nchannel; ++ch )
      {
        const double ndata = channel_counts[ch];
        if( ndata >= PEAK_FIT_MIN_CHANNEL_UNCERT )
          roi_residuals[ch] = (T(ndata) - peak_counts[ch]) / T(sqrt(ndata));
        else
          roi_residuals[ch] = peak_counts[ch]; // ad-hoc, follows PeakFitChi2Fcn
        chi2 += roi_residuals[ch] * roi_residuals[ch];

        if( irls_variances )
        {
          roi_residuals[ch] = (T(ndata) - peak_counts[ch]) / T(sqrt(irls_variances[ch]));
        }else if( m_objective == FitObjective::PoissonAllInCeres )
        {
          if( !m_fisher_variances.empty() )
            roi_residuals[ch] = (T(ndata) - peak_counts[ch]) / T(sqrt(m_fisher_variances[roi_idx][ch]));
          else
            roi_residuals[ch] = PeakFitLMObjective::poisson_deviance_residual( std::max( ndata, 0.0 ),
                                                                peak_counts[ch], deviance_floor );
        }
      }

      if( roi_models )
      {
        vector<double> &model = (*roi_models)[roi_idx];
        model.resize( nchannel );
        for( size_t ch = 0; ch < nchannel; ++ch )
        {
          if constexpr ( std::is_same_v<T, double> )
            model[ch] = peak_counts[ch];
          else
            model[ch] = peak_counts[ch].a;
        }
      }//if( roi_models )

      const double dof_val = dof_for_roi( roi_idx );
      const T chi2_dof = (dof_val > 0.0) ? (chi2 / T(dof_val)) : T(0.0);
      for( PeakType &p : peaks )
        p.set_coefficient( chi2_dof, PeakDef::Chi2DOF );

      // --- Punishment residuals (within this ROI only) ---
      // Punishment residuals immediately follow the channel residuals for this ROI.
      // Indexing: roi_residuals[nchannel + i - 1] for peaks[i] vs peaks[i-1].
      const double punishment_factor = 0.25 * static_cast<double>(nchannel) / static_cast<double>(num_roi_peaks);

      for( size_t i = 1; i < peaks.size(); ++i )
      {
        T avg_sigma;
        if constexpr ( !std::is_same_v<T, double> )
        {
          avg_sigma = T(0.5) * (peaks[i-1].sigma() + peaks[i].sigma());
        }
        else
        {
          const double s0 = peaks[i-1].gausPeak() ? peaks[i-1].sigma() : 0.25*peaks[i-1].roiWidth();
          const double s1 = peaks[i].gausPeak()   ? peaks[i].sigma()   : 0.25*peaks[i].roiWidth();
          avg_sigma = 0.5 * (s0 + s1);
        }

        if( !m_options.test( PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose ) )
        {
          const T dist = abs( peaks[i-1].mean() - peaks[i].mean() );
          T reldist = dist / avg_sigma;
          if( (reldist < 0.01) || isinf(reldist) || isnan(reldist) )
          {
            if constexpr ( !std::is_same_v<T, double> )
            {
              // Floor the *value* so the 1/reldist barrier in peaks_too_close_punishment stays
              //  finite, but keep the true derivative of reldist = |mean_i - mean_{i-1}|/sigma the
              //  division above already produced, so the fit still feels a gradient pushing the
              //  peaks apart.  If sigma->0 made reldist non-finite there is no usable derivative,
              //  so pin to a constant.
              if( isinf(reldist) || isnan(reldist) )
                reldist = T(0.01);
              else
                reldist.a = 0.01;
            }
            else
            {
              reldist = 0.01;
            }
          }

          const size_t punish_idx = nchannel + i - 1;
          assert( punish_idx < roi_residual_count( roi ) );
          // Punish peaks for being too close together; continuous + compact-support so Ceres' L-M
          //  trust region stays well-behaved (see `peaks_too_close_punishment`).
          roi_residuals[punish_idx] = peaks_too_close_punishment( reldist, punishment_factor );
        }//if( punish for peaks too close )
      }//for( size_t i = 1; i < peaks.size(); ++i )

#if( ENABLE_PUNISH_STAT_INSIG_PEAKS )
      // Punish statistically-insignificant peaks.  Lifted OUT of the proximity loop above so it runs
      //  exactly once per ROI -- it has its own peak loop; nested in the proximity loop it
      //  accumulated (npeaks-1)x via += and never ran for single-peak ROIs.  Off by default: the
      //  punishment magnitudes are heuristic and untested -- only the residual-slot layout is correct.
      if( m_options.test( PeakFitLM::PeakFitLMOptions::PunishForPeakBeingStatInsig ) )
      {
        // Use last slot of punishment residuals for the insignificance penalty
        const size_t last_punish_idx = nchannel + peaks.size() - 1;
        assert( last_punish_idx < roi_residual_count( roi ) );
        roi_residuals[last_punish_idx] = T(0.0);  // initialized once; summed below

        for( size_t pi = 0; pi < peaks.size(); ++pi )
        {
          const PeakType &pk = peaks[pi];

          if( m_options.test( PeakFitLM::PeakFitLMOptions::DoNotPunishForBeingToClose ) )
          {
            const size_t idx = nchannel + pi;
            assert( idx < roi_residual_count( roi ) );
            roi_residuals[idx] = T(0.0);
          }

          double pk_mean, pk_sigma, pk_amp;
          if constexpr ( std::is_same_v<T, double> )
          {
            pk_mean  = pk.mean();
            pk_sigma = pk.sigma();
            pk_amp   = pk.amplitude();
          }
          else
          {
            pk_mean  = pk.mean().a;
            pk_sigma = pk.sigma().a;
            pk_amp   = pk.amplitude().a;
          }

          const float lower_energy = static_cast<float>( pk_mean - 1.75*pk_sigma );
          const float upper_energy = static_cast<float>( pk_mean + 1.75*pk_sigma );
          const size_t lower_ch = std::max( roi.lower_channel, m_data->find_gamma_channel(lower_energy) );
          const size_t upper_ch = std::min( roi.upper_channel, m_data->find_gamma_channel(upper_energy) );

          double dataarea = 0.0;
          for( size_t bin = lower_ch; (bin < counts_vec.size()) && (bin <= upper_ch); ++bin )
            dataarea += counts_vec[bin];

          const double punishment_chi2 = 2.0 * static_cast<double>( 1 + upper_ch - lower_ch );
          const size_t this_punish_idx = nchannel + pi;

          if( pk_amp < std::max( 2.0*sqrt(dataarea), 1.0 ) )
          {
            if( pk_amp <= 1.0 )
            {
              if constexpr ( std::is_same_v<T, double> )
              {
                roi_residuals[this_punish_idx] += 100.0 * punishment_chi2;
              }
              else
              {
                T punish( 1.0 );
                // Normalize the Jet derivative by the amplitude value; floor the denominator so it
                //  stays finite when amp.a is ~0 (or <0) in this branch.
                const double amp_den = std::max( pk.amplitude().a, 1.0E-6 );
                punish.v = pk.amplitude().v / amp_den;
                punish *= 100.0 * punishment_chi2;
                roi_residuals[this_punish_idx] += punish;
              }
            }
            else
            {
              roi_residuals[this_punish_idx] += (T(-1.0) + 2.0*sqrt(dataarea)/pk.amplitude()) * T(punishment_chi2);
            }
          }//if( amp < 2*sqrt(dataarea) )
        }//for( size_t pi = 0; pi < peaks.size(); ++pi )
      }//if( PunishForPeakBeingStatInsig )
#endif //( ENABLE_PUNISH_STAT_INSIG_PEAKS )

      // --- Inherit user-selected options from matching starting peaks ---
      assert( roi.peaks.size() == peaks.size() );
      map<size_t,size_t> old_to_new;

      // Fixed-mean peaks keep their original index
      for( size_t idx = 0; idx < roi.peaks.size(); ++idx )
      {
        if( !roi.peaks[idx]->fitFor( PeakDef::Mean ) )
          old_to_new[idx] = idx;
      }

      for( size_t new_idx = 0; new_idx < peaks.size(); ++new_idx )
      {
        bool already = false;
        for( const auto &kv : old_to_new )
          already = (already || (kv.second == new_idx));
        if( already )
          continue;

        int closest = -1;
        const T new_mean = peaks[new_idx].mean();
        for( size_t orig_idx = 0; orig_idx < roi.peaks.size(); ++orig_idx )
        {
          if( old_to_new.count( orig_idx ) )
            continue;
          if( closest < 0 )
          {
            closest = static_cast<int>( orig_idx );
          }
          else
          {
            const double prev_m = roi.peaks[static_cast<size_t>(closest)]->mean();
            const double orig_m = roi.peaks[orig_idx]->mean();
            double nm;
            if constexpr ( std::is_same_v<T, double> )
              nm = new_mean;
            else
              nm = new_mean.a;
            if( std::abs(nm - orig_m) < std::abs(nm - prev_m) )
              closest = static_cast<int>( orig_idx );
          }
        }
        assert( closest >= 0 );
        old_to_new[static_cast<size_t>(closest)] = new_idx;
      }

      assert( old_to_new.size() == peaks.size() );
      for( const auto &mapping : old_to_new )
      {
        assert( mapping.first < peaks.size() );
        assert( mapping.second < peaks.size() );
        peaks[mapping.second].inheritUserSelectedOptions( *roi.peaks[mapping.first], true );
      }

      // VoigtPlusBortel special case: when global skew is Bortel or GaussPlusBortel
      // but an individual starting peak had VoigtPlusBortel type, remap the skew coefficients.
      // gamma_lor (SkewPar0) and R (SkewPar1, for Bortel global) are taken directly from the
      // stored original peak values (not fitted separately), since they are not in the global
      // parameter set.
      if( (m_skew_type == PeakDef::SkewType::Bortel)
          || (m_skew_type == PeakDef::SkewType::GaussPlusBortel) )
      {
        for( const auto &mapping : old_to_new )
        {
          const std::shared_ptr<const PeakDef> &orig = roi.peaks[mapping.first];
          PeakType &fitted = peaks[mapping.second];

          if( orig->skewType() != PeakDef::SkewType::VoigtPlusBortel )
            continue;

          if( m_skew_type == PeakDef::SkewType::Bortel )
          {
            // Bortel: SkewPar0=tau
            // VoigtPlusBortel: SkewPar0=gamma_lor, SkewPar1=R, SkewPar2=tau
            // roi_skew_ptr[0] is always tau for Bortel, whether shared or per-ROI.
            const T tau = roi_skew_ptr[0];
            fitted.setSkewType( PeakDef::SkewType::VoigtPlusBortel );
            fitted.set_coefficient( T(orig->coefficient(PeakDef::SkewPar0)), PeakDef::SkewPar0 ); // gamma_lor from orig
            fitted.set_coefficient( T(orig->coefficient(PeakDef::SkewPar1)), PeakDef::SkewPar1 ); // R from orig
            fitted.set_coefficient( tau, PeakDef::SkewPar2 );  // tau from skew block
            if( roi_skew_uncert_ptr )
              fitted.set_uncertainty( roi_skew_uncert_ptr[0], PeakDef::SkewPar2 );
          }
          else // GaussPlusBortel
          {
            // GaussPlusBortel: SkewPar0=R, SkewPar1=tau
            // VoigtPlusBortel: SkewPar0=gamma_lor, SkewPar1=R, SkewPar2=tau
            // roi_skew_ptr[0]/[1] are R and tau, whether shared or per-ROI.
            const T R_global   = roi_skew_ptr[0]; // GaussPlusBortel SkewPar0=R
            const T tau_global = roi_skew_ptr[1]; // GaussPlusBortel SkewPar1=tau
            fitted.setSkewType( PeakDef::SkewType::VoigtPlusBortel );
            fitted.set_coefficient( T(orig->coefficient(PeakDef::SkewPar0)), PeakDef::SkewPar0 ); // gamma_lor from orig
            fitted.set_coefficient( R_global,   PeakDef::SkewPar1 );
            fitted.set_coefficient( tau_global, PeakDef::SkewPar2 );
            if( roi_skew_uncert_ptr )
            {
              fitted.set_uncertainty( roi_skew_uncert_ptr[0], PeakDef::SkewPar1 );
              fitted.set_uncertainty( roi_skew_uncert_ptr[1], PeakDef::SkewPar2 );
            }
          }
        }//for( const auto &mapping : old_to_new )
      }//if( global skew is Bortel or GaussPlusBortel )

      return peaks;
    };//process_one_roi lambda

    // Collect all peaks from all ROIs into the final result
    vector<PeakType> all_peaks;
    all_peaks.reserve( m_total_num_peaks );

#if( PEAK_FIT_LM_PARALLEL_ROIS )
    // For >2 ROIs, evaluate each ROI on a thread from m_thread_pool.
    // The pool is sized to physical cores and lives for the lifetime of this cost function,
    // so thread creation cost is paid once at construction, not per Ceres evaluation.
    // SpecUtilsAsync::ThreadPool (GCD backend) is deliberately avoided: creating pools
    // nested inside another pool's workers causes severe delays on macOS (see PeakFit_imp.hpp).
    // Applies to both double and Jet instantiations: Jet arithmetic is ~8x more expensive
    // per channel than double, so it benefits at least as much from parallelism.
    if( nrois > 2 )
    {
      // Per-ROI residual buffers — each thread writes to its own buffer, eliminating races.
      vector<vector<T>> roi_residual_bufs( nrois );
      for( size_t i = 0; i < nrois; ++i )
        roi_residual_bufs[i].assign( roi_residual_count( m_rois[i] ), T(0.0) );

      // One promise/future pair per ROI carries the peak vector (and any exception).
      vector<std::promise<vector<PeakType>>> promises( nrois );
      vector<std::future<vector<PeakType>>>  futures;
      futures.reserve( nrois );
      for( size_t i = 0; i < nrois; ++i )
        futures.push_back( promises[i].get_future() );

      for( size_t i = 0; i < nrois; ++i )
      {
        boost::asio::post( *m_thread_pool,
          [&, i]()
          {
            try
            {
              promises[i].set_value(
                process_one_roi( i, roi_param_offsets[i], roi_residual_bufs[i].data() ) );
            }
            catch( ... )
            {
              promises[i].set_exception( std::current_exception() );
            }
          } );
      }

      // Collect results in ROI order (preserves peak ordering)
      for( size_t i = 0; i < nrois; ++i )
      {
        vector<PeakType> roi_peaks = futures[i].get(); // re-throws any stored exception
        for( PeakType &p : roi_peaks )
          all_peaks.push_back( std::move(p) );

        // Copy ROI residuals into the shared output array (sequential, no race)
        const size_t roff = roi_residual_offsets[i];
        const size_t rcount = roi_residual_count( m_rois[i] );
        for( size_t j = 0; j < rcount; ++j )
          residuals[roff + j] = roi_residual_bufs[i][j];
      }
    }else
#endif //PEAK_FIT_LM_PARALLEL_ROIS
    {
      // Serial path: used when PEAK_FIT_LM_PARALLEL_ROIS==0 or nrois<=2.
      // Each ROI writes directly into the correct slice of `residuals`.
      for( size_t roi_idx = 0; roi_idx < nrois; ++roi_idx )
      {
        T * const roi_res = residuals + roi_residual_offsets[roi_idx];
        vector<PeakType> roi_peaks = process_one_roi( roi_idx, roi_param_offsets[roi_idx], roi_res );
        for( PeakType &p : roi_peaks )
          all_peaks.push_back( std::move(p) );
      }
    }


#if( !defined(NDEBUG) && PERFORM_DEVELOPER_CHECKS && defined(CERES_PUBLIC_JET_H_) )
    const size_t nresids = number_residuals();
    for( size_t i = 0; i < nresids; ++i )
    {
      if constexpr ( std::is_same_v<T, double> )
      {
        assert( !isnan(residuals[i]) && !isinf(residuals[i]) );
      }
      else
      {
        check_jet_for_NaN( residuals[i] );
      }
    }
#endif

    return all_peaks;
  }//parametersToPeaks(...)


  /** The model counts (continuum + peaks) of every ROI's channels, at the given parameters. */
  std::vector<std::vector<double>> roi_models( const double * const params ) const
  {
    std::vector<std::vector<double>> models;
    std::vector<double> residuals( number_residuals(), 0.0 );
    parametersToPeaks<PeakDef,double>( params, nullptr, residuals.data(), nullptr, 0, &models );
    return models;
  }//roi_models(...)


  /** The statistic the IRLS passes minimize, for the given per-ROI models: the Poisson deviance of
   the ROIs being refit (`reweight`), plus the modified-Neyman chi2 of the others (which a joint fit
   with shared skew still fits by chi2).  Safeguards the IRLS passes in `run_ceres_fit(...)`.
   */
  double irls_merit( const std::vector<std::vector<double>> &models, const std::vector<bool> &reweight ) const
  {
    assert( (models.size() == m_rois.size()) && (reweight.size() == m_rois.size()) );
    const vector<float> &counts = *m_data->gamma_counts();
    double merit = 0.0;
    for( size_t r = 0; r < m_rois.size(); ++r )
    {
      const RoiInfo &roi = m_rois[r];
      for( size_t i = 0; i < models[r].size(); ++i )
      {
        const double n = static_cast<double>( counts[roi.lower_channel + i] );
        const double m = models[r][i];
        if( reweight[r] )
          merit += PeakFitLMObjective::poisson_deviance_term( std::max( 0.0, n ), std::max( m, 1.0E-12 ) );
        else
          merit += (n >= PEAK_FIT_MIN_CHANNEL_UNCERT) ? ((n - m)*(n - m)/n) : (m*m);  //as the residuals
      }
    }
    return merit;
  }//irls_merit(...)


  /** For the reweighted (IRLS) fit: sets each channel's variance to the given per-ROI model,
   floored at `sm_irls_min_variance_frac` of the ROI's mean counts per channel (but at least that
   fraction of a count), for the ROIs `reweight` marks; the others keep the chi2 weighting.
   Must not be called while the functor is being evaluated.
   */
  void set_irls_variances( const std::vector<std::vector<double>> &models, const std::vector<bool> &reweight )
  {
    assert( models.size() == m_rois.size() );
    assert( reweight.size() == m_rois.size() );

    m_irls_variances = models;
    for( size_t r = 0; r < m_rois.size(); ++r )
    {
      if( !reweight[r] )
      {
        m_irls_variances[r].clear();
        continue;
      }

      const RoiInfo &roi = m_rois[r];
      const size_t nchannel = roi.upper_channel - roi.lower_channel + 1;
      assert( m_irls_variances[r].size() == nchannel );
      const double min_variance = sm_irls_min_variance_frac
                                  * std::max( 1.0, roi.data_area / static_cast<double>(nchannel) );
      for( double &v : m_irls_variances[r] )
      {
        if( !(v >= min_variance) )  //also catches NaN
          v = min_variance;
      }
    }//for( size_t r = 0; r < m_rois.size(); ++r )
  }//set_irls_variances(...)


  /** Back to the chi2 weighting for every ROI.  Must not be called while the functor is being evaluated. */
  void clear_irls_variances()
  {
    m_irls_variances.clear();
  }


  /** Per ROI, whether its peak region is sparse enough, at the given parameters, for the Poisson
   likelihood to matter; see `sm_sparse_data_likelihood_threshold`.
   */
  std::vector<bool> sparse_rois( const double * const params ) const
  {
    std::vector<std::vector<double>> models;
    std::vector<double> residuals( number_residuals(), 0.0 );
    const std::vector<PeakDef> peaks = parametersToPeaks<PeakDef,double>( params, nullptr, residuals.data(),
                                                                           nullptr, 0, &models );
    std::vector<bool> answer( m_rois.size(), false );
    if( (models.size() != m_rois.size()) || (peaks.size() != m_total_num_peaks) )
      return answer;

    const float * const energies = m_data->channel_energies()->data();
    size_t peak_index = 0;
    for( size_t r = 0; r < m_rois.size(); ++r )
    {
      const RoiInfo &roi = m_rois[r];
      std::vector<std::pair<double,double>> mean_fwhms;
      for( size_t i = 0; i < roi.peaks.size(); ++i, ++peak_index )
        mean_fwhms.emplace_back( peaks[peak_index].mean(), peaks[peak_index].fwhm() );

      const double stat = sparse_data_statistic( energies + roi.lower_channel, models[r].data(),
                                                 models[r].size(), mean_fwhms );
      answer[r] = (stat > sm_sparse_data_likelihood_threshold);
    }

    return answer;
  }//sparse_rois(...)


  /** Adds, to the uncertainty of every quantity the linear sub-solve provides - the peak amplitudes,
   and the polynomial continuum coefficients - the part due to the non-linear parameters' own
   uncertainty, so the reported uncertainties are marginal rather than conditional on the fitted
   means, widths and skew.

   With x the linear and theta the non-linear parameters, the block inverse of the full
   Gauss-Newton covariance gives
     Cov(x) = Cov(x | theta) + G Cov(theta) G^T,   G = dx/dtheta,
   where Cov(x | theta) is what the linear solve reports, and Cov(theta) is what Ceres computes from
   the residuals with x solved out (variable projection).  G is the derivative of the linear solution with respect to theta,
   obtained by evaluating with Jets (seeded eight parameters at a time).

   The propagation is linear, so it is only as good as Cov(theta).  Ceres' covariance knows nothing
   of the parameter bounds, and for a degenerate fit (e.g., an insignificant peak whose width sits at
   a bound, or a step coefficient with nothing to step) it can be enormous; each parameter's standard
   deviation is therefore capped at that of a uniform distribution over its allowed range,
   (upper - lower)/sqrt(12), keeping the correlations.

   `peaks` must be the `parametersToPeaks<PeakDef,double>(...)` output for `params`; amplitudes that
   are Ceres parameters (a ROI off the LLS path) already have marginal uncertainties and are left
   alone, as are the peak-CDF step coefficients.
   */
  void add_nonlinear_uncertainty_to_linear_pars( const double * const params,
                                                 std::vector<double> row_major_covariance,
                                                 const ProblemSetup &prob_setup,
                                                 std::vector<PeakDef> &peaks ) const
  {
    using Jet8 = ceres::Jet<double,8>;
    const size_t npar = m_num_parameters;
    if( row_major_covariance.size() != npar*npar )
      return;

    // Only quantities the linear sub-solve provides need this.
    if( std::none_of( begin(m_rois), end(m_rois), []( const RoiInfo &r ){ return r.use_lls_for_cont; } ) )
      return;

    std::vector<size_t> free_pars;
    for( size_t k = 0; k < npar; ++k )
    {
      if( row_major_covariance[k*npar + k] > 0.0 )
        free_pars.push_back( k );
    }
    if( free_pars.empty() || (peaks.size() != m_total_num_peaks) )
      return;

    // Cap each standard deviation at the spread of a uniform distribution over the bounds.
    for( const size_t k : free_pars )
    {
      if( (k >= prob_setup.m_lower_bounds.size()) || !prob_setup.m_lower_bounds[k].has_value()
         || !prob_setup.m_upper_bounds[k].has_value() )
        continue;
      const double max_sd = (*prob_setup.m_upper_bounds[k] - *prob_setup.m_lower_bounds[k]) / std::sqrt( 12.0 );
      const double sd = std::sqrt( row_major_covariance[k*npar + k] );
      if( (max_sd > 0.0) && (sd > max_sd) )
      {
        const double factor = max_sd / sd;
        for( size_t j = 0; j < npar; ++j )
        {
          row_major_covariance[k*npar + j] *= factor;
          row_major_covariance[j*npar + k] *= factor;
        }
      }
    }//for( const size_t k : free_pars )

    const size_t nfree = free_pars.size();
    const size_t max_cont = 5;  // PeakContinuumImp holds up to five parameters
    std::vector<std::vector<double>> amp_grad( peaks.size(), std::vector<double>( nfree, 0.0 ) );
    std::vector<std::vector<std::vector<double>>> cont_grad( peaks.size(),
                     std::vector<std::vector<double>>( max_cont, std::vector<double>( nfree, 0.0 ) ) );

    std::vector<Jet8> jet_params( npar );
    std::vector<Jet8> jet_residuals( number_residuals() );
    for( size_t first = 0; first < nfree; first += 8 )
    {
      for( size_t k = 0; k < npar; ++k )
        jet_params[k] = Jet8( params[k] );
      for( size_t j = first; (j < nfree) && (j < first + 8); ++j )
        jet_params[free_pars[j]].v[j - first] = 1.0;

      const std::vector<PeakDefImpWithCont<Jet8>> jet_peaks
          = parametersToPeaks<PeakDefImpWithCont<Jet8>,Jet8>( jet_params.data(), nullptr, jet_residuals.data() );
      if( jet_peaks.size() != peaks.size() )
        return;

      for( size_t i = 0; i < jet_peaks.size(); ++i )
      {
        const PeakDefImpWithCont<Jet8> &jp = jet_peaks[i];
        for( size_t j = first; (j < nfree) && (j < first + 8); ++j )
        {
          amp_grad[i][j] = jp.amplitude().v[j - first];
          if( jp.m_continuum )
          {
            for( size_t c = 0; c < max_cont; ++c )
              cont_grad[i][c][j] = jp.m_continuum->parameters()[c].v[j - first];
          }
        }
      }//for( peaks )
    }//for( blocks of eight free parameters )

    // g^T Cov(theta) g, over the free parameters.
    const auto propagated_variance = [&]( const std::vector<double> &g ) -> double {
      double var = 0.0;
      for( size_t a = 0; a < nfree; ++a )
      {
        if( g[a] == 0.0 )
          continue;
        for( size_t b = 0; b < nfree; ++b )
          var += g[a] * row_major_covariance[free_pars[a]*npar + free_pars[b]] * g[b];
      }
      return (var > 0.0) ? var : 0.0;
    };

    size_t peak_index = 0;
    for( const RoiInfo &roi : m_rois )
    {
      const bool amps_from_lls = (roi_amp_parameter_count( roi ) == 0);
      const size_t num_lls_cont = roi.use_lls_for_cont ? PeakContinuum::num_linear_fit_pars( roi.offset_type ) : 0;

      for( size_t n = 0; n < roi.peaks.size(); ++n, ++peak_index )
      {
        PeakDef &peak = peaks[peak_index];
        if( amps_from_lls )
        {
          const double extra = propagated_variance( amp_grad[peak_index] );
          if( extra > 0.0 )
          {
            const double cond = std::max( peak.amplitudeUncert(), 0.0 );
            peak.setAmplitudeUncert( std::sqrt( cond*cond + extra ) );
          }
        }

        // The ROI's peaks share one continuum; update it once, from its first peak.
        if( (n == 0) && num_lls_cont )
        {
          const std::shared_ptr<PeakContinuum> cont = peak.getContinuum();
          const std::vector<double> uncerts = cont->uncertainties();
          for( size_t c = 0; (c < num_lls_cont) && (c < max_cont) && (c < uncerts.size()); ++c )
          {
            const double extra = propagated_variance( cont_grad[peak_index][c] );
            if( extra > 0.0 )
              cont->setPolynomialUncert( c, std::sqrt( uncerts[c]*uncerts[c] + extra ) );
          }
        }
      }//for( peaks in this ROI )
    }//for( ROIs )
  }//add_nonlinear_uncertainty_to_linear_pars(...)


  /** Sets (or, with an empty argument, clears) `m_fisher_variances` from per-ROI models, floored as
   the deviance residuals floor the model.  Must not be called while the functor is being evaluated.
   */
  void set_fisher_variances( const std::vector<std::vector<double>> &models )
  {
    m_fisher_variances = models;
    for( size_t r = 0; r < m_fisher_variances.size(); ++r )
    {
      const RoiInfo &roi = m_rois[r];
      const size_t nchannel = roi.upper_channel - roi.lower_channel + 1;
      const double floor = sm_deviance_min_expected_frac
                           * std::max( 1.0, roi.data_area / static_cast<double>(nchannel) );
      for( double &v : m_fisher_variances[r] )
      {
        if( !(v >= floor) )
          v = floor;
      }
    }
  }//set_fisher_variances(...)


  template<typename T>
  bool operator()( T const *const *parameters, T *residuals ) const
  {
    m_ncalls += 1;

    try
    {
      T const * const pars = parameters[0];
      const T * const uncertainties = nullptr;
      if constexpr ( std::is_same_v<T, double> )
      {
        parametersToPeaks<PeakDef,double>( pars, uncertainties, residuals );
      }else
      {
        parametersToPeaks<PeakDefImpWithCont<T>,T>( pars, uncertainties, residuals );
      }
    }catch( std::exception &e )
    {
      cerr << "PeakFitDiffCostFunction: caught exception: '" << e.what() << "'" << endl;
      return false;
    }

    return true;
  }//operator()


  static PeakDef::SkewType skew_type_from_prev_peaks( const vector<shared_ptr<const PeakDef> > &inpeaks )
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
  }//void skew_type_from_prev_peaks(...)


  /** Set up the shared skew parameter block (params[0..skew_block_size-1]).
   Skew params are at the front of the global parameter vector.
   Layout:
     params[0..num_skew-1]:        base (lower-anchor) values for all skew params
     params[num_skew..num_skew+M-1]: upper-anchor values for energy-dep params only
       (only present when m_fit_skew_energy_dependence == true)
   When IndependentSkewValues is set, this function is a no-op (skew is per-ROI).
  */
  void setup_skew_parameters( double *pars,
                              vector<int> &constant_parameters,
                              vector<std::optional<double>> &lower_bounds,
                              vector<std::optional<double>> &upper_bounds,
                              const vector<shared_ptr<const PeakDef>> &inpeaks ) const
  {
    // With IndependentSkewValues, there is no shared skew block; each ROI sets up its own.
    if( m_options.test( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
      return;

    const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skew_type );

    if( num_skew_pars == 0 )
      return;

    vector<bool> fit_parameter( num_skew_pars, true );
    vector<double> starting_value( num_skew_pars, 0.0 );
    vector<double> lower_values( num_skew_pars, 0.0 );
    vector<double> upper_values( num_skew_pars, 0.0 );

    for( size_t i = 0; i < num_skew_pars; ++i )
    {
      const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
      double lower, upper, start, dx;
      const bool use = PeakDef::skew_parameter_range( m_skew_type, ct, lower, upper, start, dx );
      assert( use );
      if( !use )
        throw logic_error( "Inconsistent skew par val" );
      starting_value[i] = start;
      lower_values[i]   = lower;
      upper_values[i]   = upper;
    }

    // Use values from the first matching peak (matching global skew type)
    for( const auto &p : inpeaks )
    {
      if( p->skewType() != m_skew_type )
        continue;

      for( size_t i = 0; i < num_skew_pars; ++i )
      {
        const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
        fit_parameter[i] = p->fitFor( ct );
        double val = p->coefficient( ct );

        if( IsInf(val) || IsNan(val) || (val < lower_values[i]) || (val > upper_values[i]) )
          val = starting_value[i];

        // Sanity-clamp Crystal Ball power-law params that can drift high
        switch( m_skew_type )
        {
          case PeakDef::NumSkewType:
            assert( 0 );
            // Fall through to NoSkew for non-debug builds
          case PeakDef::NoSkew:   case PeakDef::Bortel:
          case PeakDef::DoubleBortel: case PeakDef::GaussPlusBortel:
          case PeakDef::GaussExp: case PeakDef::ExpGaussExp:
          case PeakDef::GadrasGeneric: case PeakDef::GadrasCZT:
            break;

          case PeakDef::CrystalBall:
          case PeakDef::DoubleSidedCrystalBall:
          {
            switch( ct )
            {
              case PeakDef::Mean:           case PeakDef::Sigma:
              case PeakDef::GaussAmplitude: case PeakDef::NumCoefficientTypes:
              case PeakDef::Chi2DOF:
              case PeakDef::SkewPar4:       case PeakDef::SkewPar5:
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
            }
            break;
          }

          case PeakDef::VoigtPlusBortel:
            break;
        }//switch( m_skew_type )

        starting_value[i] = val;
      }
    }//for( const auto &p : inpeaks )

    // For small/medium refinement, restrict skew range near starting value to prevent large shifts;
    // otherwise allow the full parameter range so fitting from default starting values can succeed.
    const bool restrict_skew_range = m_options.test( PeakFitLM::PeakFitLMOptions::SmallAmplitudeRefinementOnly )
                                     || m_options.test( PeakFitLM::PeakFitLMOptions::MediumAmplitudeRefinementOnly );

    // Set up the base (lower-anchor) skew parameters at params[0..num_skew_pars-1]
    for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
    {
      pars[skew_index] = starting_value[skew_index];

      if( fit_parameter[skew_index] )
      {
        if( restrict_skew_range )
        {
          // A window centred on the starting value, floored at a quarter of the parameter's full
          //  range so a zero starting value doesnt collapse it to a point - see
          //  `LinearProblemSubSolveChi2Fcn::addSkewParameters`.
          //
          // Expressed as an OFFSET from the starting value rather than a multiple of it: scaling by
          //  the value silently inverts for a negative parameter (lower bound above upper bound,
          //  i.e. an infeasible starting point handed to Ceres).  No skew parameter is negative
          //  today, so this is latent - but it bit immediately when a log-scaled skew coordinate
          //  was trialled, and cost an afternoon to find.
          const double par_range = upper_values[skew_index] - lower_values[skew_index];
          const double half_window = (std::max)( 0.125*par_range,
                                                0.5*fabs(starting_value[skew_index]) );
          lower_bounds[skew_index] = (std::max)( lower_values[skew_index],
                                                starting_value[skew_index] - half_window );
          upper_bounds[skew_index] = (std::min)( upper_values[skew_index],
                                                starting_value[skew_index] + half_window );
        }else
        {
          lower_bounds[skew_index] = lower_values[skew_index];
          upper_bounds[skew_index] = upper_values[skew_index];
        }

        assert( lower_bounds[skew_index] < upper_bounds[skew_index] );
      }
      else
      {
        constant_parameters.push_back( static_cast<int>(skew_index) );
      }
    }

    // For energy-dependent multi-ROI fitting, set up upper-anchor params
    if( m_fit_skew_energy_dependence )
    {
      size_t upper_idx = num_skew_pars;
      for( size_t i = 0; i < num_skew_pars; ++i )
      {
        const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
        if( !PeakDef::is_energy_dependent( m_skew_type, ct ) )
          continue;

        // Initialize upper-anchor to same value as lower-anchor
        pars[upper_idx] = starting_value[i];
        if( fit_parameter[i] )
        {
          lower_bounds[upper_idx] = lower_bounds[i];
          upper_bounds[upper_idx] = upper_bounds[i];
        }
        else
        {
          constant_parameters.push_back( static_cast<int>(upper_idx) );
        }
        upper_idx += 1;
      }
    }//if( m_fit_skew_energy_dependence )
  }//setup_skew_parameters(...)


  /** Bound on the dimensionless peak-CDF step coefficients.

   The Ceres variable is the step's continuum-density change divided by the continuum level, so
   +-2 says the step may not move the continuum by more than twice its own level.  A fit that
   wants more than that has gone wrong, not found a steeper step.
   */
  static constexpr double sm_cdf_step_bound = 2.0;


  /** Pre-solves this ROI's peak-CDF step coefficients by least-squares, and overwrites `seeds`
   (in dimensionless units) if it succeeds.

   With every peak amplitude held at its input value the step coefficients are linear, so
   `PeakFit::fit_continuum(...)` can solve them directly.  Failure - a singular design matrix, say
   - is a legitimate "the data says nothing about the step here" outcome, so we quietly keep the
   incoming seeds.
   */
  void warm_start_cdf_step( const RoiInfo &roi, vector<double> &seeds ) const
  {
    if( !m_data || !m_data->num_gamma_channels() || (roi.upper_channel <= roi.lower_channel) )
      return;

    try
    {
      const shared_ptr<const vector<float>> channel_energies = m_data->channel_energies();
      const shared_ptr<const vector<float>> gamma_counts = m_data->gamma_counts();
      if( !channel_energies || !gamma_counts
         || ((roi.upper_channel + 1) >= channel_energies->size())
         || (roi.upper_channel >= gamma_counts->size()) )
        return;

      // Note the "+ 1": the basis must sit on exactly the channel grid the residual uses, or the
      //  seed is anchored differently from the fit.
      const size_t nchannel = 1 + roi.upper_channel - roi.lower_channel;
      const float * const energies = &((*channel_energies)[roi.lower_channel]);
      const float * const counts = &((*gamma_counts)[roi.lower_channel]);

      vector<PeakDef> fixed_amp_peaks;
      fixed_amp_peaks.reserve( roi.peaks.size() );
      for( const shared_ptr<const PeakDef> &p : roi.peaks )
        fixed_amp_peaks.push_back( *p );

      const size_t num_cont_pars = PeakContinuum::num_parameters( roi.offset_type );
      const size_t num_poly = PeakContinuum::num_linear_fit_pars( roi.offset_type );
      vector<double> cont_coeffs( num_cont_pars, 0.0 );
      vector<double> dummy_peak_counts( nchannel, 0.0 );

      // `roi.ref_energy`, not the continuum's own reference energy - make_rois may have moved it,
      //  and a step paired with a polynomial about a different origin is not the same step.
      PeakFit::fit_continuum( energies, counts, static_cast<const float *>(nullptr),
                              nchannel, roi.offset_type, roi.ref_energy,
                              fixed_amp_peaks, false,
                              cont_coeffs.data(), dummy_peak_counts.data() );

      for( size_t k = 0; (k < seeds.size()) && ((num_poly + k) < cont_coeffs.size()); ++k )
      {
        const double solved = cont_coeffs[num_poly + k] / cdf_step_par_scale( roi, k );
        if( std::isfinite(solved) )
          seeds[k] = solved;
      }
    }catch( std::exception & )
    {
      // Keep the incoming seeds.
    }
  }//void warm_start_cdf_step( const RoiInfo &roi, vector<double> &seeds ) const


  /** Set up parameters for a single ROI at param_offset in the global parameter array. */
  void setup_roi_parameters( const RoiInfo &roi,
                             const size_t param_offset,
                             double *pars,
                             vector<int> &constant_parameters,
                             vector<std::optional<double>> &lower_bounds,
                             vector<std::optional<double>> &upper_bounds ) const
  {
    const size_t num_fit_cont = roi_cont_parameter_count( roi );
    const size_t num_sigmas_fit = roi_sigma_parameter_count( roi );
    const size_t num_amps_fit = roi_amp_parameter_count( roi );
    const double range = roi.upper_energy - roi.lower_energy;

    // Peak-CDF step coefficients occupy the first slots of this ROI's parameter block, held in
    //  dimensionless form (see RoiInfo::cdf_step_scale).
    if( roi.use_lls_for_cont && num_fit_cont )
    {
      const shared_ptr<const PeakContinuum> initial_continuum = roi.peaks.front()->continuum();
      const vector<double> &cont_pars = initial_continuum->parameters();
      const size_t num_poly = PeakContinuum::num_linear_fit_pars( roi.offset_type );

      // make_rois may have moved the reference energy, and BiLinearStepCDF's step is a function of
      //  E' - so translate the stored coefficients into the ROI's frame before using them:
      //    s0 + s1*(E - old_ref) == (s0 + s1*(new_ref - old_ref)) + s1*(E - new_ref)
      vector<double> stored_pars( num_fit_cont, 0.0 );
      for( size_t k = 0; k < num_fit_cont; ++k )
        stored_pars[k] = ((num_poly + k) < cont_pars.size()) ? cont_pars[num_poly + k] : 0.0;

      const double d_ref = roi.ref_energy - initial_continuum->referenceEnergy();
      if( (num_fit_cont > 1) && std::isfinite(d_ref) && (d_ref != 0.0) )
        stored_pars[0] += stored_pars[1] * d_ref;

      vector<double> seeds( num_fit_cont, 0.0 ), stored_seeds( num_fit_cont, 0.0 );
      bool need_warm_start = false;
      for( size_t k = 0; k < num_fit_cont; ++k )
      {
        const double stored = stored_pars[k];
        seeds[k] = stored / cdf_step_par_scale( roi, k );
        stored_seeds[k] = seeds[k];

        // A stored value of exactly zero means nobody has fit this step yet - e.g. the user just
        //  switched continuum type.  chi2 is shallow in this parameter, so starting from zero
        //  usually *ends* near zero, turning a FlatStepCDF into a plain Constant.
        need_warm_start |= (!std::isfinite(seeds[k]) || (seeds[k] == 0.0)
                            || (seeds[k] < -sm_cdf_step_bound) || (seeds[k] > sm_cdf_step_bound));
      }

      if( need_warm_start )
        warm_start_cdf_step( roi, seeds );

      const vector<bool> cont_fit_for = initial_continuum->fitForParameter();

      for( size_t k = 0; k < num_fit_cont; ++k )
      {
        // Ceres rejects an infeasible starting point outright, so clamp rather than widen - the
        //  whole value of the bound is that it is physically meaningful.
        pars[param_offset + k] = std::max( -sm_cdf_step_bound,
                                           std::min( sm_cdf_step_bound, seeds[k] ) );

        // A step coefficient the caller pinned stays where it started; the polynomial terms are
        //  still LLS-solved, so a pinned step does not cost us the LLS path.
        if( ((num_poly + k) < cont_fit_for.size()) && !cont_fit_for[num_poly + k] )
        {
          // The caller's value, not whatever the warm start solved for.
          pars[param_offset + k] = std::isfinite(stored_seeds[k]) ? stored_seeds[k] : 0.0;
          constant_parameters.push_back( static_cast<int>(param_offset + k) );
          continue;
        }

        lower_bounds[param_offset + k] = -sm_cdf_step_bound;
        upper_bounds[param_offset + k] = sm_cdf_step_bound;
      }
    }
    // Continuum parameters (if not LLS)
    else if( !roi.use_lls_for_cont )
    {
      const shared_ptr<const PeakContinuum> initial_continuum = roi.peaks.front()->continuum();
      assert( initial_continuum );
      const vector<double> &cont_pars    = initial_continuum->parameters();
      const vector<bool>    par_fit_for  = initial_continuum->fitForParameter();
      assert( cont_pars.size() == num_fit_cont );

      for( size_t i = 0; i < num_fit_cont; ++i )
      {
        pars[param_offset + i] = cont_pars[i];
        if( !par_fit_for[i] )
          constant_parameters.push_back( static_cast<int>(param_offset + i) );
      }

      // Bound the step coefficients here too - on this path they are raw Ceres parameters handed
      //  straight to `setParameters(...)`, so they stay in keV^-1 rather than being rescaled.
      //  Same +-2 meaning as the LLS path, just expressed in physical units.
      const size_t num_poly = PeakContinuum::num_linear_fit_pars( roi.offset_type );
      for( size_t k = 0; (num_poly + k) < num_fit_cont; ++k )
      {
        const size_t index = param_offset + num_poly + k;
        if( !par_fit_for[num_poly + k] )
          continue;   // a pinned parameter is constant; bounding it would only risk infeasibility

        const double limit = sm_cdf_step_bound * cdf_step_par_scale( roi, k );
        pars[index] = std::max( -limit, std::min( limit, pars[index] ) );
        lower_bounds[index] = -limit;
        upper_bounds[index] = limit;
      }
    }

    // Compute mean range bin width for sigma constraints
    const size_t midbin = m_data->find_gamma_channel( static_cast<float>(0.5*(roi.lower_energy + roi.upper_energy)) );
    const float binwidth = m_data->gamma_channel_width( midbin );
    const double avrg_bin_width = [&](){
      const double left_energy  = m_data->gamma_channel_lower( roi.lower_channel );
      const double right_energy = m_data->gamma_channel_upper( roi.upper_channel );
      return (right_energy - left_energy) / static_cast<double>( 1 + roi.upper_channel - roi.lower_channel );
    }();

    double minsigma = binwidth, max_input_sigma = binwidth, min_input_sigma = binwidth;
    double maxsigma = 0.5*range;
    // The input widths themselves - the two above are floored at a channel width, which made the
    //  Small/Medium FWHM refinement's lower bound (relative to the narrowest input) a fraction of a
    //  channel for any peak wider than one, so a width could shrink without limit (a NaI line at
    //  2 MeV refit to a third of the resolution).  HPGe keeps the historical bounds until the change
    //  is measured and accepted there (on Detective-X it moved 52 of 192 fits, one moderate line lost).
    const bool width_bounds_from_input = (m_det_type != PeakFitUtils::CoarseResolutionType::High);
    double narrowest_input_sigma = std::numeric_limits<double>::infinity(), widest_input_sigma = 0.0;

    // Mean parameters and compute minsigma/maxsigma
    for( size_t i = 0; i < roi.peaks.size(); ++i )
    {
      const std::shared_ptr<const PeakDef> &peak = roi.peaks[i];
      const double mean  = peak->mean();
      const double sigma = peak->sigma();
      const size_t mean_par_index = param_offset + num_fit_cont + num_sigmas_fit + i;

      const double rel_mean = 0.5 + (mean - roi.lower_energy) / range;
      pars[mean_par_index] = rel_mean;

      if( !peak->fitFor( PeakDef::Mean ) )
      {
        constant_parameters.push_back( static_cast<int>(mean_par_index) );
      }
      else
      {
        lower_bounds[mean_par_index] = 0.5;
        upper_bounds[mean_par_index] = 1.5;

        if( m_options.test( PeakFitLM::PeakFitLMOptions::MediumAmplitudeRefinementOnly ) )
        {
          lower_bounds[mean_par_index] = 0.5 + (mean - 0.5*sigma - roi.lower_energy) / range;
          upper_bounds[mean_par_index] = 0.5 + (mean + 0.5*sigma - roi.lower_energy) / range;
          assert( rel_mean >= *lower_bounds[mean_par_index] );
          assert( rel_mean <= *upper_bounds[mean_par_index] );
        }

        if( m_options.test( PeakFitLM::PeakFitLMOptions::SmallAmplitudeRefinementOnly ) )
        {
          lower_bounds[mean_par_index] = 0.5 + (mean - 0.15*sigma - roi.lower_energy) / range;
          upper_bounds[mean_par_index] = 0.5 + (mean + 0.15*sigma - roi.lower_energy) / range;
          assert( rel_mean >= *lower_bounds[mean_par_index] );
          assert( rel_mean <= *upper_bounds[mean_par_index] );
        }

        if( rel_mean <= *lower_bounds[mean_par_index] )
          lower_bounds[mean_par_index] = rel_mean - std::max( 0.1, fabs(0.25*rel_mean) );
        if( rel_mean >= *upper_bounds[mean_par_index] )
          upper_bounds[mean_par_index] = 1.2*rel_mean;

        assert( rel_mean >= *lower_bounds[mean_par_index] );
        assert( rel_mean <= *upper_bounds[mean_par_index] );
      }

      min_input_sigma = std::min( min_input_sigma, sigma );
      max_input_sigma = std::max( max_input_sigma, sigma );
      narrowest_input_sigma = std::min( narrowest_input_sigma, sigma );
      widest_input_sigma = std::max( widest_input_sigma, sigma );

      if( !peak->fitFor( PeakDef::Sigma ) )
      {
        minsigma = std::min( minsigma, sigma );
        maxsigma = std::max( maxsigma, sigma );
      }
      else
      {
        float lowersigma, uppersigma;
        expected_peak_width_limits( mean, m_det_type, m_data, lowersigma, uppersigma );
        if( i == 0 )
          minsigma = lowersigma;
        if( i == (roi.peaks.size() - 1) )
          maxsigma = uppersigma;
      }
    }//for each peak

    // Peak amplitude parameters; only present when the LLS is not solving them (see
    //  roi_amp_parameter_count).  Scaled by the input amplitude so the variable starts at 1.
    if( num_amps_fit )
    {
      size_t fit_amp_num = 0;
      for( size_t i = 0; i < roi.peaks.size(); ++i )
      {
        if( !roi.peaks[i]->fitFor( PeakDef::GaussAmplitude ) )
          continue;

        const size_t amp_index = param_offset + num_fit_cont + num_sigmas_fit
                                 + roi.peaks.size() + fit_amp_num;
        const double amp_scale = roi_amp_par_scale( roi, i );

        // The all-in-Ceres objective lets a peak's area go negative, as the LLS does for every other
        //  objective, so a peak consistent with zero is not clipped (which biases it up); the
        //  deviance residuals stay finite where the model dips below zero.  Otherwise a peak area
        //  cannot be negative.  The bound's magnitude just keeps a runaway in check.
        const bool allow_negative = (m_objective == FitObjective::PoissonAllInCeres);
        pars[amp_index] = roi.peaks[i]->amplitude() / amp_scale;
        if( !std::isfinite(pars[amp_index]) || (!allow_negative && (pars[amp_index] < 0.0)) )
          pars[amp_index] = 0.0;

        upper_bounds[amp_index] = std::max( 100.0, 10.0*fabs(pars[amp_index]) );
        lower_bounds[amp_index] = allow_negative ? -upper_bounds[amp_index].value() : 0.0;

        fit_amp_num += 1;
      }//for( size_t i = 0; i < roi.peaks.size(); ++i )
    }//if( num_amps_fit )

    // Sigma parameters
    if( num_sigmas_fit == 0 )
    {
      // No sigma parameters to set up
    }
    else if( m_options.test( PeakFitLM::PeakFitLMOptions::AllPeakFwhmIndependent ) )
    {
      // One parameter per peak with fitted sigma
      size_t fit_sigma_index = 0;
      for( const std::shared_ptr<const PeakDef> &peak : roi.peaks )
      {
        if( !peak->fitFor( PeakDef::Sigma ) )
          continue;

        const double sigma = peak->sigma();
        const size_t sigma_par_idx = param_offset + num_fit_cont + fit_sigma_index;

        pars[sigma_par_idx] = sigma / roi.max_initial_sigma;
        lower_bounds[sigma_par_idx] = (0.5*std::min(minsigma,min_input_sigma)) / roi.max_initial_sigma;
        upper_bounds[sigma_par_idx] = (1.5*std::max(maxsigma,max_input_sigma)) / roi.max_initial_sigma;

        if( m_options.test( PeakFitLM::PeakFitLMOptions::MediumFwhmRefinementOnly ) )
        {
          lower_bounds[sigma_par_idx] = 0.5 * sigma / roi.max_initial_sigma;
          upper_bounds[sigma_par_idx] = 1.5 * sigma / roi.max_initial_sigma;
        }
        if( m_options.test( PeakFitLM::PeakFitLMOptions::SmallFwhmRefinementOnly ) )
        {
          lower_bounds[sigma_par_idx] = 0.85 * sigma / roi.max_initial_sigma;
          upper_bounds[sigma_par_idx] = 1.15 * sigma / roi.max_initial_sigma;
        }

        assert( pars[sigma_par_idx] >= *lower_bounds[sigma_par_idx] );
        assert( pars[sigma_par_idx] <= *upper_bounds[sigma_par_idx] );
        fit_sigma_index += 1;
      }
      assert( fit_sigma_index == num_sigmas_fit );
    }
    else
    {
      // One or two shared sigma parameters.
      // The parameter represents sigma / roi.max_initial_sigma.
      // roi.max_initial_sigma has a 1.0 keV floor, so the starting value may be
      // less than 1.0 when the actual peak sigmas are smaller than 1.0 keV.
      const size_t fit_sigma_idx = param_offset + num_fit_cont;
      pars[fit_sigma_idx] = max_input_sigma / roi.max_initial_sigma;

      lower_bounds[fit_sigma_idx] = (0.5*std::min(minsigma, 0.75*min_input_sigma)) / roi.max_initial_sigma;
      upper_bounds[fit_sigma_idx] = (1.5*std::max(maxsigma, 1.25*max_input_sigma)) / roi.max_initial_sigma;

      assert( pars[fit_sigma_idx] >= *lower_bounds[fit_sigma_idx] );
      assert( pars[fit_sigma_idx] <= *upper_bounds[fit_sigma_idx] );

      const double refine_min_sigma = width_bounds_from_input ? narrowest_input_sigma : min_input_sigma;
      const double refine_max_sigma = width_bounds_from_input ? widest_input_sigma : max_input_sigma;
      if( m_options.test( PeakFitLM::PeakFitLMOptions::MediumFwhmRefinementOnly ) )
      {
        lower_bounds[fit_sigma_idx] = (0.5*refine_min_sigma) / roi.max_initial_sigma;
        upper_bounds[fit_sigma_idx] = (1.5*refine_max_sigma) / roi.max_initial_sigma;
        if( width_bounds_from_input )
          pars[fit_sigma_idx] = std::min( pars[fit_sigma_idx], *upper_bounds[fit_sigma_idx] );
        assert( pars[fit_sigma_idx] >= *lower_bounds[fit_sigma_idx] );
        assert( pars[fit_sigma_idx] <= *upper_bounds[fit_sigma_idx] );
      }
      if( m_options.test( PeakFitLM::PeakFitLMOptions::SmallFwhmRefinementOnly ) )
      {
        lower_bounds[fit_sigma_idx] = (0.85*refine_min_sigma) / roi.max_initial_sigma;
        upper_bounds[fit_sigma_idx] = (1.15*refine_max_sigma) / roi.max_initial_sigma;
        if( width_bounds_from_input )
          pars[fit_sigma_idx] = std::min( pars[fit_sigma_idx], *upper_bounds[fit_sigma_idx] );
        assert( pars[fit_sigma_idx] >= *lower_bounds[fit_sigma_idx] );
        assert( pars[fit_sigma_idx] <= *upper_bounds[fit_sigma_idx] );
      }

      // Enforce minimum of 0.525 channels per sigma (narrowest commercially observed detector)
      bool any_peak_fixed_narrow = false, all_widths_fixed = true;
      for( const auto &peak : roi.peaks )
      {
        const bool fit_width = peak->fitFor( PeakDef::Sigma );
        any_peak_fixed_narrow |= (!fit_width && (peak->sigma() < 0.525*avrg_bin_width));
        all_widths_fixed = (all_widths_fixed && !fit_width);
      }

      if( !any_peak_fixed_narrow )
      {
        const double rel_bin_width = avrg_bin_width / roi.max_initial_sigma;
        lower_bounds[fit_sigma_idx] = std::max( 0.525*rel_bin_width, lower_bounds[fit_sigma_idx].value() );
        upper_bounds[fit_sigma_idx] = std::max( 0.525*rel_bin_width, upper_bounds[fit_sigma_idx].value() );
        pars[fit_sigma_idx]         = std::max( 0.525*rel_bin_width, pars[fit_sigma_idx] );
      }

      assert( pars[fit_sigma_idx] >= *lower_bounds[fit_sigma_idx] );
      assert( pars[fit_sigma_idx] <= *upper_bounds[fit_sigma_idx] );
      assert( all_widths_fixed || (lower_bounds[fit_sigma_idx].has_value()
              && (lower_bounds[fit_sigma_idx].value() >= 0.525*(avrg_bin_width/roi.max_initial_sigma))) );

      if( num_sigmas_fit > 1 )
      {
        // Second sigma param = multiplier for the high-energy end; range +-20%
        const size_t upper_sigma_idx = fit_sigma_idx + 1;
        pars[upper_sigma_idx] = 1.0;
        lower_bounds[upper_sigma_idx] = 0.8;
        upper_bounds[upper_sigma_idx] = 1.2;
        assert( pars[upper_sigma_idx] >= *lower_bounds[upper_sigma_idx] );
        assert( pars[upper_sigma_idx] <= *upper_bounds[upper_sigma_idx] );
      }
    }//sigma parameter setup

    // Per-ROI skew parameters (only when IndependentSkewValues option is set)
    if( m_options.test( PeakFitLM::PeakFitLMOptions::IndependentSkewValues ) )
    {
      const size_t num_skew_pars = PeakDef::num_skew_parameters( m_skew_type );
      if( num_skew_pars > 0 )
      {
        // Determine default ranges and starting values for each skew parameter
        vector<bool>   fit_parameter( num_skew_pars, false ); // OR'd across matching peaks below
        vector<double> starting_value( num_skew_pars, 0.0 );
        vector<double> lower_values( num_skew_pars, 0.0 );
        vector<double> upper_values( num_skew_pars, 0.0 );

        for( size_t i = 0; i < num_skew_pars; ++i )
        {
          const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
          double lower, upper, start, dx;
          const bool use = PeakDef::skew_parameter_range( m_skew_type, ct, lower, upper, start, dx );
          assert( use );
          if( !use )
            throw logic_error( "Inconsistent skew par val (IndependentSkewValues)" );
          starting_value[i] = start;
          lower_values[i]   = lower;
          upper_values[i]   = upper;
        }

        // A skew parameter is fitted if ANY matching-type peak in the ROI has fitFor==true for it;
        // it is held constant only when ALL matching peaks agree it should be fixed.
        // fit_parameter[] starts false; we OR in each matching peak's fitFor flag.
        // Initial coefficient values come from the first matching peak.
        bool found_first = false;
        for( const auto &p : roi.peaks )
        {
          if( p->skewType() != m_skew_type )
            continue;

          for( size_t i = 0; i < num_skew_pars; ++i )
          {
            const auto ct = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );

            // OR across peaks: mark as fit if any peak requests it
            if( p->fitFor( ct ) )
              fit_parameter[i] = true;

            if( !found_first )
            {
              // Coefficient starting value: use first matching peak only
              double val = p->coefficient( ct );

              if( IsInf(val) || IsNan(val) || (val < lower_values[i]) || (val > upper_values[i]) )
                val = starting_value[i];

              // Sanity-clamp Crystal Ball power-law params that can drift high
              switch( m_skew_type )
              {
                case PeakDef::NumSkewType:
                  assert( 0 );
                  // Fall through to NoSkew for non-debug builds
                case PeakDef::NoSkew:   case PeakDef::Bortel:
                case PeakDef::DoubleBortel: case PeakDef::GaussPlusBortel:
                case PeakDef::GaussExp: case PeakDef::ExpGaussExp:
                case PeakDef::GadrasGeneric: case PeakDef::GadrasCZT:
                  break;

                case PeakDef::CrystalBall:
                case PeakDef::DoubleSidedCrystalBall:
                {
                  switch( ct )
                  {
                    case PeakDef::Mean:           case PeakDef::Sigma:
                    case PeakDef::GaussAmplitude: case PeakDef::NumCoefficientTypes:
                    case PeakDef::Chi2DOF:
                    case PeakDef::SkewPar4:       case PeakDef::SkewPar5:
                    case PeakDef::SkewPar0:
                    case PeakDef::SkewPar2:
                      if( val > 3.0 )
                        val = starting_value[i];
                      break;
                    case PeakDef::SkewPar1:
                    case PeakDef::SkewPar3:
                      if( val > 6.0 )
                        val = starting_value[i];
                      break;
                  }
                  break;
                }

                case PeakDef::VoigtPlusBortel:
                  break;
              }//switch( m_skew_type )

              starting_value[i] = val;
            }//if( !found_first )
          }//for( size_t i = 0; i < num_skew_pars; ++i )

          found_first = true;
        }//for( const auto &p : roi.peaks )

        // Write the per-ROI skew params into the global parameter array
        // Must match process_one_roi's `roi_skew_ptr`: the skew block sits after the amplitude
        //  block, which is only present when the LLS is not solving the amplitudes.
        const size_t skew_base_idx = param_offset + num_fit_cont + num_sigmas_fit
                                     + roi.peaks.size() + num_amps_fit;
        const bool restrict_skew_range = m_options.test( PeakFitLM::PeakFitLMOptions::SmallAmplitudeRefinementOnly )
                                         || m_options.test( PeakFitLM::PeakFitLMOptions::MediumAmplitudeRefinementOnly );
        for( size_t skew_index = 0; skew_index < num_skew_pars; ++skew_index )
        {
          const size_t abs_idx = skew_base_idx + skew_index;
          pars[abs_idx] = starting_value[skew_index];

          if( fit_parameter[skew_index] )
          {
            if( restrict_skew_range )
            {
              // Floor the window at a quarter of the parameters full range, so a zero starting
              //  value doesnt collapse it to a point.
              const double par_range = upper_values[skew_index] - lower_values[skew_index];
              lower_bounds[abs_idx] = std::max( lower_values[skew_index], 0.5*starting_value[skew_index] );
              upper_bounds[abs_idx] = std::min( upper_values[skew_index],
                                    std::max( 1.5*fabs(starting_value[skew_index]), 0.25*par_range ) );
            }else
            {
              lower_bounds[abs_idx] = lower_values[skew_index];
              upper_bounds[abs_idx] = upper_values[skew_index];
            }

            assert( lower_bounds[abs_idx] < upper_bounds[abs_idx] );
          }
          else
          {
            constant_parameters.push_back( static_cast<int>(abs_idx) );
          }
        }
      }//if( num_skew_pars > 0 )
    }//if( IndependentSkewValues )
  }//setup_roi_parameters(...)


  ProblemSetup get_problem_setup() const
  {
    const size_t num_fit_pars = number_parameters();

    vector<int> constant_parameters;
    vector<double> parameters( num_fit_pars, 0.0 );
    vector<std::optional<double>> lower_bounds( num_fit_pars );
    vector<std::optional<double>> upper_bounds( num_fit_pars );

    // Shared skew parameters at the front
    setup_skew_parameters( &parameters[0], constant_parameters, lower_bounds, upper_bounds,
                           all_starting_peaks() );

    // Per-ROI parameters following the skew block
    size_t param_offset = skew_parameter_count();
    for( const RoiInfo &roi : m_rois )
    {
      setup_roi_parameters( roi, param_offset, &parameters[0], constant_parameters,
                            lower_bounds, upper_bounds );
      param_offset += roi_parameter_count( roi );
    }

    assert( param_offset == num_fit_pars );

    ProblemSetup prob_setup;
    prob_setup.m_parameters         = parameters;
    prob_setup.m_constant_parameters = constant_parameters;
    prob_setup.m_lower_bounds        = lower_bounds;
    prob_setup.m_upper_bounds        = upper_bounds;

    return prob_setup;
  }//get_problem_setup()

public:
  // NOTE: declaration order must match initialization dependency order.
  //  m_rois depends on m_data, m_objective; m_external_continuum depends on m_rois;
  //  m_fit_skew_energy_dependence depends on m_rois; m_num_parameters depends on
  //  m_fit_skew_energy_dependence; m_thread_pool depends on m_rois.
  const std::shared_ptr<const SpecUtils::Measurement> m_data;
  const PeakDef::SkewType m_skew_type;
  const PeakFitUtils::CoarseResolutionType m_det_type;
  const Wt::WFlags<PeakFitLMOptions> m_options;
  const FitObjective m_objective;

  mutable std::atomic<unsigned int> m_ncalls;

  /** For a ROI being refit by IRLS: per channel, the variance each channel's residual (and the
   linear solve) is weighted by - the model of the previous fit pass, floored.  Empty for a ROI fit by
   chi2, and for the first pass, which is the chi2 fit.  Only changed between solves
   (`set_irls_variances`), never during an evaluation, so concurrent evaluations of the ROIs stay safe.
   */
  std::vector<std::vector<double>> m_irls_variances;

  /** For `FitObjective::PoissonAllInCeres`: when not empty (per ROI, per channel), the channel
   residuals are `(data - model)/sqrt(variance)` instead of the signed root deviance - set to the
   fitted model while Ceres computes the covariance, so `J^T J` is the expected Fisher information
   (the root-deviance residuals' `J^T J` gives empty channels half their information).  Like
   `m_irls_variances`, only changed between solves.
   */
  std::vector<std::vector<double>> m_fisher_variances;

  const std::vector<RoiInfo> m_rois;
  const size_t m_total_num_peaks;
  const bool m_fit_skew_energy_dependence;
  const double m_skew_anchor_lower_energy;
  const double m_skew_anchor_upper_energy;
  const std::shared_ptr<const SpecUtils::Measurement> m_external_continuum;
  const size_t m_num_parameters;
  const size_t m_num_residuals;

#if( PEAK_FIT_LM_PARALLEL_ROIS )
  // Non-null only when m_rois.size() > 2; sized to min(nrois, physical_cores).
  const std::unique_ptr<boost::asio::thread_pool> m_thread_pool;
#endif
};//struct PeakFitDiffCostFunction


/** Results of solving a PeakFitDiffCostFunction: the fitted peaks, and the raw fit
 parameter/uncertainty/covariance data.
 */
struct CeresFitResult
{
  vector<PeakDef> final_peaks;
  vector<double> parameters;
  vector<double> uncertainties;
  vector<double> row_major_covariance;
  size_t num_fit_pars = 0;

  /** Some ROI was fit by Poisson likelihood (refit by IRLS, or the all-in-Ceres objective). */
  bool by_likelihood = false;

  /** The solve was aborted by the cancel flag; nothing else is valid. */
  bool cancelled = false;
};


static CeresFitResult run_ceres_fit( PeakFitDiffCostFunction &cost_functor, const size_t total_num_peaks,
                                     const std::shared_ptr<const std::atomic<bool>> &cancel_flag = nullptr );


/** Evaluates the model at the starting parameters, for a problem where every Ceres parameter is held
 constant - e.g., a Peak Editor refit with centroid and FWHM fixed, where only the amplitudes and
 continuum vary, and the linear least-squares in `parametersToPeaks(...)` solves those directly.
 Ceres would have nothing to minimize, and rejects the zero initial trust-region radius such a
 problem gives it.  Zero uncertainties/covariance are what Ceres reports for constant parameters, and
 they make `parametersToPeaks(...)` keep the input peaks' mean/FWHM uncertainties.
 */
static CeresFitResult evaluate_without_free_parameters( const PeakFitDiffCostFunction &cost_functor,
                                         const PeakFitDiffCostFunction::ProblemSetup &prob_setup )
{
  const size_t num_fit_pars = prob_setup.m_parameters.size();
  assert( prob_setup.m_constant_parameters.size() == num_fit_pars );

  CeresFitResult result;
  result.num_fit_pars = num_fit_pars;
  result.parameters = prob_setup.m_parameters;
  result.uncertainties.resize( num_fit_pars, 0.0 );
  result.row_major_covariance.resize( num_fit_pars * num_fit_pars, 0.0 );

  vector<double> residuals( cost_functor.number_residuals(), 0.0 );
  result.final_peaks = cost_functor.parametersToPeaks<PeakDef,double>( result.parameters.data(),
                         result.uncertainties.data(), residuals.data(),
                         result.row_major_covariance.data(), num_fit_pars );

  return result;
}//evaluate_without_free_parameters(...)


/** `fit_peaks_in_roi_LM(...)`, also setting `by_likelihood` to whether the ROI was fit by Poisson
 likelihood.
 */
static vector<shared_ptr<const PeakDef>> fit_peaks_in_roi_imp( const vector<shared_ptr<const PeakDef>> &coFitPeaks,
                                                      const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                                      const PeakFitUtils::CoarseResolutionType det_type,
                                                      const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options,
                                                      const std::shared_ptr<const std::atomic<bool>> &cancel_flag,
                                                      bool &by_likelihood )
{
  by_likelihood = false;

  /** For this first go, we will have Ceres fit for things.

   This is only tested to the "seems to work, for simple cases" level.

   In the future:
     - Need to deal with errors during solving, and during Covariance computation
     - Deal with setting number of threads reasonably.
     - Try out using a ceres::LossFunction
     - Optimize/pcik the Ceres-based paramaters to get reasonable fits.
   */

  if( coFitPeaks.empty() )
    throw runtime_error( "fit_peaks_in_roi_LM: empty input peaks." );

  // Early bail-out: if cancellation has already been signaled before we set up
  //  the Ceres problem, don't waste time building it just to abort on the
  //  first iteration.
  if( cancel_flag && cancel_flag->load( std::memory_order_relaxed ) )
    return {};

  // Lets make sure all coFitPeaks share a continuum.
  for( size_t i = 1; i < coFitPeaks.size(); ++i )
  {
    if( coFitPeaks[i]->continuum() != coFitPeaks[0]->continuum() )
      throw runtime_error( "fit_peak_for_user_click_LM: input peaks all must share a single continuum" );
  }//for( size_t i = 1; i < coFitPeaks.size(); ++i )

  try
  {
    // TODO: check that repeaded calls to this function wont cause adding/removing channels due to find_gamma_channel(...) rounding type things - e.g., do we need to subtract off a tiny bit from roiUpperEnergy

    const shared_ptr<const PeakContinuum> input_continuum = coFitPeaks[0]->continuum();

    const float roiLowerEnergy = static_cast<float>( input_continuum->lowerEnergy() );
    const float roiUpperEnergy = static_cast<float>( input_continuum->upperEnergy() );
    assert( roiLowerEnergy < roiUpperEnergy );
    if( roiLowerEnergy >= roiUpperEnergy )
      throw runtime_error( "Invalid energy range (" + std::to_string(roiLowerEnergy) + ", " + std::to_string(roiUpperEnergy) + ")" );


    const size_t lower_channel = dataH->find_gamma_channel(roiLowerEnergy);
    const size_t upper_channel = dataH->find_gamma_channel(roiUpperEnergy);
    const double roi_lower = dataH->gamma_channel_lower( lower_channel );
    const double roi_upper = dataH->gamma_channel_upper( upper_channel );
    
    assert( (lower_channel < 2) || (roi_lower <= roiLowerEnergy) );
    assert( (roi_upper >= roiUpperEnergy) || ((upper_channel + 3) >= dataH->num_gamma_channels()) );

    if( coFitPeaks.size() >= (upper_channel - lower_channel) ) //PeakFitDiffCostFunction will check for this in an assert
      throw runtime_error( "Invalid energy channels (" + std::to_string(lower_channel) + ", "
                          + std::to_string(upper_channel) + ") for " + std::to_string(coFitPeaks.size()) + " peaks." );
    assert( coFitPeaks.size() < (upper_channel - lower_channel) );

    const PeakContinuum::OffsetType offset_type = input_continuum->type();
    const size_t num_continuum_pars = PeakContinuum::num_parameters( offset_type );

    const double prev_ref_energy = input_continuum->referenceEnergy();
    const bool prev_ref_energy_valid = ((prev_ref_energy >= roiLowerEnergy) && (prev_ref_energy <= roiUpperEnergy));
    const double reference_energy = prev_ref_energy_valid ? prev_ref_energy : roiLowerEnergy;

    const vector<bool> parFitFors = input_continuum->fitForParameter();
    assert( parFitFors.size() <= num_continuum_pars );

    const PeakDef::SkewType skew_type = PeakFitDiffCostFunction::skew_type_from_prev_peaks( coFitPeaks );

    // With `SPARSE_DATA_LIKELIHOOD_USE_CERES`, the likelihood fit is all-in-Ceres, which starts from
    //  the chi2 fit (its amplitudes and continuum are Ceres parameters, so need good starting values);
    //  by default, that chi2 fit also decides whether the ROI is sparse enough to be refit.
    const Wt::WFlags<PeakFitLMOptions> eff_options = effective_options( fit_options, dataH );
    const bool forced = eff_options.test( PeakFitLMOptions::ForcePoissonLikelihood );
#if( SPARSE_DATA_LIKELIHOOD_USE_CERES )
    const bool likelihood_by_ceres = !eff_options.test( PeakFitLMOptions::NoSparseDataLikelihood )
                                     && (offset_type != PeakContinuum::OffsetType::External);
#else
    const bool likelihood_by_ceres = false;
#endif
    FitObjective objective = FitObjective::NeymanChi2;
    vector<shared_ptr<const PeakDef>> start_peaks = coFitPeaks;
    if( likelihood_by_ceres )
    {
      bool chi2_by_likelihood = false;
      const vector<shared_ptr<const PeakDef>> chi2_peaks = fit_peaks_in_roi_imp( coFitPeaks, dataH, det_type,
                                                chi2_only_options( fit_options ), cancel_flag, chi2_by_likelihood );
      if( chi2_peaks.size() != coFitPeaks.size() )
        return chi2_peaks;  //cancelled

      if( !forced )
      {
        if( sparse_data_statistic( chi2_peaks, dataH ) <= sm_sparse_data_likelihood_threshold )
          return chi2_peaks;
        tl_fit_objective_diagnostics.sparse_rois += 1;
      }

      objective = FitObjective::PoissonAllInCeres;
      start_peaks = chi2_peaks;
    }//if( likelihood_by_ceres )

    // The solve, its retries, any reweighting (IRLS) passes, and the covariance are shared with
    //  `fit_peaks_in_spectrum_LM(...)`.  A failed IRLS refit keeps the chi2 fit (see
    //  `run_ceres_fit(...)`), so a failure here is the chi2 solve itself - except for all-in-Ceres,
    //  which then falls back to the chi2 fit it started from.
    CeresFitResult fit;
    try
    {
      PeakFitDiffCostFunction cost_functor( dataH, start_peaks, roiLowerEnergy, roiUpperEnergy,
                                            reference_energy, skew_type, det_type, fit_options, objective );
      fit = run_ceres_fit( cost_functor, coFitPeaks.size(), cancel_flag );
    }catch( std::exception & )
    {
      if( objective == FitObjective::NeymanChi2 )
        throw;

      tl_fit_objective_diagnostics.fell_back_to_chi2 = true;
      return start_peaks;
    }//try / catch

    if( fit.cancelled )
      return {};

    vector<shared_ptr<const PeakDef>> results( fit.final_peaks.size() );
    for( size_t i = 0; i < fit.final_peaks.size(); ++i )
      results[i] = make_shared<PeakDef>( std::move( fit.final_peaks[i] ) );
    by_likelihood = fit.by_likelihood;

    return results;
  }catch( std::exception &e )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    cout << "fit_peak_for_user_click_LM caught: " << e.what() << endl;
#endif
    //assert( 0 );
    throw;
  }//try / catch

  assert( 0 );
  throw runtime_error( "shouldnt have gotten here" );
  return {};
}//fit_peaks_in_roi_imp(...)


/** All peaks passed in must share a PeakContinuum.
 */
vector<shared_ptr<const PeakDef>> fit_peaks_in_roi_LM( const vector<shared_ptr<const PeakDef>> coFitPeaks,
                                                      const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                                                      const PeakFitUtils::CoarseResolutionType det_type,
                                                      const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options,
                                                      std::shared_ptr<const std::atomic<bool>> cancel_flag )
{
  bool by_likelihood = false;
  return fit_peaks_in_roi_imp( coFitPeaks, dataH, det_type, fit_options, cancel_flag, by_likelihood );
}//fit_peaks_in_roi_LM(...)



void fit_peak_for_user_click_LM( PeakShrdVec &results,
                             const std::shared_ptr<const SpecUtils::Measurement> &dataH,
                             const vector<shared_ptr<const PeakDef>> &coFitPeaksInput,
                             const double mean0, const double sigma0,
                             const double area0,
                             const float roiLowerEnergy,
                             const float roiUpperEnergy,
                             const std::shared_ptr<const PeakFitDetPrefs> &fitPrefs,
                             const std::shared_ptr<const DetectorPeakResponse> &drf,
                             std::shared_ptr<const std::atomic<bool>> cancel_flag )
{
  vector<shared_ptr<const PeakDef>> coFitPeaks = coFitPeaksInput;

  const PeakDef::SkewType skew_type = PeakFitDiffCostFunction::skew_type_from_prev_peaks( coFitPeaks );
  const size_t num_skew = PeakDef::num_skew_parameters( skew_type );

  std::shared_ptr<PeakDef> candidatepeak = std::make_shared<PeakDef>(mean0, sigma0, area0);

  // If the load-time prefs left m_det_type == Unknown, fall back to classifying from the
  // spectrum so the fitter has a usable det_type (otherwise narrow HPGe peaks get fit with
  // low-res defaults and rejected).
  PeakFitUtils::CoarseResolutionType det_type
    = PeakFitUtils::effective_det_type( fitPrefs, dataH, nullptr );

  // Apply skew type from fitPrefs to the candidate peak, overriding any
  //  skew type inferred from existing nearby peaks.
  assert( fitPrefs );
  if( fitPrefs )
  {
    candidatepeak->setSkewType( fitPrefs->m_peak_skew_type );

    const size_t num_prefs_skew
      = PeakDef::num_skew_parameters( fitPrefs->m_peak_skew_type );
    for( size_t i = 0; i < num_prefs_skew; ++i )
    {
      const PeakDef::CoefficientType ct
        = PeakDef::CoefficientType( PeakDef::CoefficientType::SkewPar0 + i );

      if( fitPrefs->m_lower_energy_skew[i].has_value() )
      {
        double val = fitPrefs->m_lower_energy_skew[i].value();

        if( PeakDef::is_energy_dependent( fitPrefs->m_peak_skew_type, ct )
           && fitPrefs->m_upper_energy_skew[i].has_value()
           && dataH && dataH->num_gamma_channels() > 0 )
        {
          const double lower_energy = dataH->gamma_channel_lower( 0 );
          const double upper_energy
            = dataH->gamma_channel_upper( dataH->num_gamma_channels() - 1 );
          const double frac = (mean0 - lower_energy) / (upper_energy - lower_energy);
          val += frac * (fitPrefs->m_upper_energy_skew[i].value() - val);
        }

        candidatepeak->set_coefficient( val, ct );
        candidatepeak->setFitFor( ct, false );
      }//if( skew param value specified )
    }//for( loop over skew parameters )

    // Apply FWHM method from preferences
    if( (fitPrefs->m_fwhm_method != PeakFitDetPrefs::FwhmMethod::Normal)
       && drf && drf->hasResolutionInfo() )
    {
      const double drf_sigma
        = drf->peakResolutionSigma( static_cast<float>( mean0 ) );
      candidatepeak->setSigma( drf_sigma );

      if( fitPrefs->m_fwhm_method == PeakFitDetPrefs::FwhmMethod::DetFwhm )
        candidatepeak->setFitFor( PeakDef::CoefficientType::Sigma, false );
    }//if( apply FWHM from DRF )
  }else
  {
    const PeakFitSpec::SpecClassType type = PeakFitSpec::initial_lowres_highres_classify( dataH );

    switch( type )
    {
      case PeakFitSpec::SpecClassType::LowOrMedRes:
        det_type = PeakFitUtils::CoarseResolutionType::LowOrMedRes;
      break;

      case PeakFitSpec::SpecClassType::High:
        det_type = PeakFitUtils::CoarseResolutionType::High;
      break;

      case PeakFitSpec::SpecClassType::Unknown:
        det_type = PeakFitUtils::CoarseResolutionType::Unknown;
      break;
    }
  }//if( fitPrefs )

  const bool isHPGe = (det_type == PeakFitUtils::CoarseResolutionType::High);

  //The below should probably go off the number of bins in the ROI
  const size_t num_fit_peaks = coFitPeaksInput.size() + 1;
  PeakContinuum::OffsetType offset_type = (num_fit_peaks < (isHPGe ? 3 : 2)) ? PeakContinuum::Linear : PeakContinuum::Quadratic;

  // TODO: remove isHPGe usage above and use det_type directly

  if( coFitPeaks.size() )
    offset_type = std::max( offset_type, coFitPeaks[0]->continuum()->type() );


  // Make sure all initial peaks share a continuum
  if( coFitPeaks.size() )
  {
    shared_ptr<PeakContinuum> initial_continuum = make_shared<PeakContinuum>( *coFitPeaks[0]->continuum() );
    initial_continuum->setRange( roiLowerEnergy, roiUpperEnergy );
    candidatepeak->setContinuum( initial_continuum );
    for( size_t i = 0; i < coFitPeaks.size(); ++i )
    {
      auto newpeak = make_shared<PeakDef>( *coFitPeaks[i] );
      newpeak->setContinuum( initial_continuum );
      coFitPeaks[i] = newpeak;
    }

    if( initial_continuum->type() != offset_type )
    {
      const vector<bool> parFitFors = initial_continuum->fitForParameter();
      assert( parFitFors.size() <= PeakContinuum::num_parameters( initial_continuum->type() ) );

      bool any_continuum_par_fixed = false;
      for( size_t fit_for_index = 0; fit_for_index < parFitFors.size(); ++fit_for_index )
        any_continuum_par_fixed |= !parFitFors[fit_for_index];

      if( !any_continuum_par_fixed )
      {
        const double prev_ref_energy = initial_continuum->referenceEnergy();
        const bool prev_ref_energy_valid = ((prev_ref_energy >= roiLowerEnergy) && (prev_ref_energy <= roiUpperEnergy));
        const double reference_energy = prev_ref_energy_valid ? prev_ref_energy : roiLowerEnergy;

        initial_continuum->calc_linear_continuum_eqn( dataH, reference_energy, roiLowerEnergy, roiUpperEnergy, 2, 2 );
      }
      initial_continuum->setType( offset_type );
    }
  }else
  {
    shared_ptr<PeakContinuum> initial_continuum = candidatepeak->continuum();
    initial_continuum->setRange( roiLowerEnergy, roiUpperEnergy );
    const double reference_energy = roiLowerEnergy;
    initial_continuum->calc_linear_continuum_eqn( dataH, reference_energy, roiLowerEnergy, roiUpperEnergy, 2, 2 );
    initial_continuum->setType( offset_type );
  }//if( coFitPeaks.size() )

  coFitPeaks.push_back( candidatepeak );
  std::sort( coFitPeaks.begin(), coFitPeaks.end(), &PeakDef::lessThanByMeanShrdPtr );


  // The peak-search acceptance tests this fit feeds (e.g., the area-significance cut in
  //  `check_highres_single_peak_fit(...)`) were tuned on chi2 fits with area uncertainties
  //  conditional on the fit mean and width, so keep those here until the tests are retuned.
  Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options
                    = Wt::WFlags<PeakFitLM::PeakFitLMOptions>( PeakFitLM::ConditionalAreaUncertainties )
                      | PeakFitLM::NoSparseDataLikelihood;
  if( fitPrefs
     && (fitPrefs->m_fwhm_method == PeakFitDetPrefs::FwhmMethod::DetPlusRefine)
     && drf && drf->hasResolutionInfo() )
  {
    fit_options |= PeakFitLM::SmallFwhmRefinementOnly;
  }

  try
  {
    results = fit_peaks_in_roi_LM( coFitPeaks, dataH, det_type, fit_options, cancel_flag );
  }catch( std::exception &e )
  {
    results.clear();
  }
}//void fit_peak_for_user_click_LM(...)


void fit_peaks_LM( vector<shared_ptr<const PeakDef>> &results,
                  const vector<shared_ptr<const PeakDef>> input_peaks,
                  shared_ptr<const SpecUtils::Measurement> data,
                  const double stat_threshold,
                  const double hypothesis_threshold,
                  const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options,
                  const PeakFitUtils::CoarseResolutionType det_type ) throw()
{
  // Relax the significance test when the caller is refining an existing fit.  `test` is a
  //  bitwise AND, and the Medium/Small composites cover both their amplitude and FWHM bits, so
  //  these two tests catch all four refinement-only options.
  const bool is_refit = fit_options.test( PeakFitLM::PeakFitLMOptions::MediumRefinementOnly )
                        || fit_options.test( PeakFitLM::PeakFitLMOptions::SmallRefinementOnly );

  try
  {
    // Check all input peaks share a ROI
    for( size_t i = 1; i < input_peaks.size(); ++i )
    {
      if( input_peaks[i]->continuum() != input_peaks[0]->continuum() )
        throw runtime_error( "fit_peaks_LM: all input peaks must share ROI." );
    }

    results.clear();

    //We have to separate out non-gaussian peaks since they cant enter the
    //  fitting methods
    vector<shared_ptr<const PeakDef>> near_peaks, datadefined_peaks;

    shared_ptr<const SpecUtils::Measurement> ext_continuum;

    near_peaks.reserve( input_peaks.size() );

    for( const shared_ptr<const PeakDef> &p : input_peaks )
    {
      if( p->gausPeak() )
        near_peaks.push_back( p );
      else
        datadefined_peaks.push_back( p );

      if( !!p->continuum()->externalContinuum() )
        ext_continuum = p->continuum()->externalContinuum();
    }//for( const PeakDef &p : all_near_peaks )


    if( near_peaks.empty() )
      return;

    // Completely replace contents of `near_peaks` with new peaks, and new `PeakContinuum`
    //  so we dont affect the input peaks.
    local_unique_copy_continuum( near_peaks );

    //Need to make sure near_peaks and fixedpeaks are all gaussian (if not
    //  seperate them out, and add them in later).  If fitpeaks is non-gaussian
    //  ignore it or throw an exception.

    double lowx( 0.0 ), highx( 0.0 );

    {
      const shared_ptr<const PeakDef> &lowgaus = near_peaks.front();
      const shared_ptr<const PeakDef> &highgaus = near_peaks.back();

      double dummy = 0.0;
      const bool isHPGe_for_roi = (det_type == PeakFitUtils::CoarseResolutionType::High);
      findROIEnergyLimits( lowx, dummy, *lowgaus, data, isHPGe_for_roi );
      findROIEnergyLimits( dummy, highx, *highgaus, data, isHPGe_for_roi );
    }


    results = fit_peaks_in_roi_LM( near_peaks, data, det_type, fit_options );

    // I'm pretty sure things have been sorted, but lets check
    assert( std::is_sorted(begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr ) );

    //
    //set_chi2_dof( data, results, 0, near_peaks.size() );


    for( size_t i = 1; i <= results.size(); ++i ) //Note weird convntion of index
    {
      const shared_ptr<const PeakDef> &peak = (results[i-1]);

      // Dont remove peaks whos amplitudes we arent fitting
      if( !peak->fitFor(PeakDef::GaussAmplitude) )
        continue;

      // Dont enforce a significance test if we are refitting the peak, and
      //  Sigma and Mean are fixed - the user probably knows what they
      //  are are doing.
      if( is_refit
         && !peak->fitFor(PeakDef::Sigma)
         && !peak->fitFor(PeakDef::Mean) )
      {
        continue;
      }

      const double num_sigma = (peak->amplitude() / peak->amplitudeUncert());
      const bool is_sig = (stat_threshold <= 0.0) || (num_sigma >= stat_threshold);

      const double dummy_stat_thresh = 0.0;
      const bool significant = chi2_significance_test( *peak, dummy_stat_thresh, hypothesis_threshold, {}, data );
      if( !is_sig || !significant )
      {
#if( PRINT_DEBUG_INFO_FOR_PEAK_SEARCH_FIT_LEVEL > 0 )
        DebugLog(cerr) << "\tPeak at mean=" << peak->mean()
        << "is being discarded for not being significant"
        << "\n";
#endif
        results.erase( results.begin() + --i );
      }//if( !significant )
    }//for( size_t i = 1; i < fitpeaks.size(); ++i )

    bool removed_peak = false;
    for( size_t i = 1; i < results.size(); ++i ) //Note weird convention of index
    {
      const shared_ptr<const PeakDef> &this_peak = results[i-1];
      const shared_ptr<const PeakDef> &next_peak = results[i+1-1];

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
          results.erase( begin(results) + i - 1 );
        else
          results.erase( begin(results) + i );

        i = i - 1; //incase we have multiple close peaks in a row
      }//if( (mean_diff/min_sigma) < 1.0 ) / else
    }//for( size_t i = 1; i < fitpeaks.size(); ++i )

    if( removed_peak )
      results = fit_peaks_in_roi_LM( results, data, det_type, fit_options );

    if( datadefined_peaks.size() )
    {
      results.insert( end(results), begin(datadefined_peaks), end(datadefined_peaks) );
      std::sort( begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr );
    }

    return;
  }catch( std::exception &e )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO || !defined(NDEBUG) )
    cerr << "fit_peaks_LM: caught exception '" << e.what() << "'" << endl;
#endif
    results.clear();
    return;
  }//try / catch


  //We will only reach here if there was no exception, so since never expect
  //  this to actually happen, just assign the results to be same as the input
  assert( 0 );
  results = input_peaks;
}//vector<PeakDef> fit_peaks_LM(...);


vector<shared_ptr<const PeakDef>> fit_peaks_in_range_LM( const double x0, const double x1,
                                      const double ncausalitysigma,
                                      const double stat_threshold,
                                      const double hypothesis_threshold,
                                      const std::vector<std::shared_ptr<const PeakDef>> input_peaks,
                                      const std::shared_ptr<const SpecUtils::Measurement> data,
                                      const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options,
                                      const PeakFitUtils::CoarseResolutionType det_type )
{
  if( !data || (x1 < x0) )
    return input_peaks;

  vector<shared_ptr<const PeakDef>> all_peaks = input_peaks;
  std::sort( begin(all_peaks), end(all_peaks), &PeakDef::lessThanByMeanShrdPtr );

  const vector<shared_ptr<const PeakDef>> peaks_in_range = peaksInRange( x0, x1, ncausalitysigma, all_peaks );
  vector<shared_ptr<const PeakDef>> peaks_not_in_range = all_peaks;

  peaks_not_in_range.erase( std::remove_if( begin(peaks_not_in_range), end(peaks_not_in_range),
    [&peaks_in_range]( const shared_ptr<const PeakDef> &p ) -> bool {
      // Return false if the peak is in the energy range of interest
      return std::find(begin(peaks_in_range), end(peaks_in_range), p) != end(peaks_in_range);
    } ),
    end(peaks_not_in_range) );

  assert( (peaks_not_in_range.size() + peaks_in_range.size()) == input_peaks.size() );


  // The returned pointers all point to the original peaks.
  const vector<vector<shared_ptr<const PeakDef>>> seperated_peaks
                                          = causilyDisconnectedPeaks( ncausalitysigma, false, peaks_in_range );

  assert( ([seperated_peaks](){ size_t i = 0; for(auto &p:seperated_peaks) i += p.size(); return i; })() == peaks_in_range.size() );

  //Fit each of the ranges
  vector<vector<shared_ptr<const PeakDef>>> fit_peak_ranges( seperated_peaks.size() );
  SpecUtilsAsync::ThreadPool threadpool;
  for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  {
    threadpool.post( [&fit_peak_ranges, &seperated_peaks, data, stat_threshold,
                      hypothesis_threshold, fit_options, det_type, peakn](){
      fit_peaks_LM( fit_peak_ranges[peakn], seperated_peaks[peakn],
                    data, stat_threshold, hypothesis_threshold, fit_options, det_type );
    } );
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )
  threadpool.join();

  vector<shared_ptr<const PeakDef>> results = peaks_not_in_range;

  //put the fit peaks back into 'input_peaks' so we can return all the peaks
  //  passed in, not just ones in X range of interest
  for( size_t peakn = 0; peakn < fit_peak_ranges.size(); ++peakn )
  {
    const vector<shared_ptr<const PeakDef>> &fit_peaks = fit_peak_ranges[peakn];
    results.insert( end(results), begin(fit_peaks), end(fit_peaks) );
  }//for( size_t peakn = 0; peakn < fit_peak_ranges.size(); ++peakn )

  std::sort( begin(results), end(results), &PeakDef::lessThanByMeanShrdPtr );

  //Now make sure peaks from two previously causally disconnected regions
  //  didn't migrate towards each other, causing the regions to become
  //  causally connected now
  bool migration = false;
  for( size_t peakn = 1; peakn < fit_peak_ranges.size(); ++peakn )
  {
    if( fit_peak_ranges[peakn-1].empty() || fit_peak_ranges[peakn].empty() )
      continue;

    const shared_ptr<const PeakDef> &last_peak = fit_peak_ranges[peakn-1].back();
    const shared_ptr<const PeakDef> &this_peak = fit_peak_ranges[peakn][0];
    migration |= PeakDef::causilyConnected( *last_peak, *this_peak, ncausalitysigma, false );
  }//for( size_t peakn = 0; peakn < seperated_peaks.size(); ++peakn )

  if( migration )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    cerr << "fitPeaksInRange(...)\n\tWarning: Migration happened!" << endl;
#endif

    return fit_peaks_in_range_LM( x0, x1, ncausalitysigma,
                           stat_threshold, hypothesis_threshold,
                                 results, data, fit_options, det_type );
  }//if( migration )

  //  cout << "Fit took: " << timer.format() << endl;

  return results;
}//vector<shared_ptr<const PeakDef>> fit_peaks_in_range_LM(...)

/** The statistic a fit of `peaks` (which must share one continuum) minimized, over channels
 [lower_channel, upper_channel]: the Poisson deviance if the ROI was fit by likelihood, otherwise the
 chi2 of `chi2_for_region(...)`.
 */
static double roi_fit_statistic( const vector<shared_ptr<const PeakDef>> &peaks,
                                 const shared_ptr<const SpecUtils::Measurement> &data,
                                 const int lower_channel, const int upper_channel,
                                 const bool by_likelihood )
{
  if( !by_likelihood || peaks.empty() || (upper_channel < lower_channel) )
    return chi2_for_region( peaks, data, lower_channel, upper_channel );

  const size_t ch0 = static_cast<size_t>( lower_channel );
  const size_t nchan = static_cast<size_t>( upper_channel - lower_channel ) + 1;
  const vector<float> &counts = *data->gamma_counts();
  const vector<double> model = roi_model_counts( peaks, data, ch0, nchan );

  double deviance = 0.0;
  for( size_t i = 0; i < nchan; ++i )
  {
    const double n = std::max( 0.0, static_cast<double>( counts[ch0 + i] ) );
    deviance += PeakFitLMObjective::poisson_deviance_term( n, std::max( model[i], 1.0E-12 ) );
  }

  return deviance;
}//roi_fit_statistic(...)


std::vector<std::shared_ptr<const PeakDef>> refitPeaksThatShareROI_LM(
                                   const std::shared_ptr<const SpecUtils::Measurement> &data,
                                   const std::shared_ptr<const DetectorPeakResponse> &detector,
                                   const std::vector<std::shared_ptr<const PeakDef>> &inpeaks,
                                   const PeakFitUtils::CoarseResolutionType det_type,
                                   const Wt::WFlags<PeakFitLM::PeakFitLMOptions> fit_options )
{
  vector<shared_ptr<const PeakDef>> answer;

  try
  {
    if( inpeaks.empty() )
      return answer;

    std::shared_ptr<const PeakContinuum> origCont = inpeaks[0]->continuum();

    for( const shared_ptr<const PeakDef> &p : inpeaks )
      if( origCont != p->continuum() )
        throw runtime_error( "refitPeaksThatShareROI_LM: all input peaks must share a ROI" );

    bool by_likelihood = false;
    answer = fit_peaks_in_roi_imp( inpeaks, data, det_type, fit_options, nullptr, by_likelihood );


    //now we need to go through and make sure the peaks we're adding are both
    //  significant, and improve the chi2/dof from before the fit.
    const double lx = origCont->lowerEnergy();
    const double ux = origCont->upperEnergy();
    const int lower_channel = static_cast<int>( data->find_gamma_channel( lx ) );
    const int upper_channel = static_cast<int>( data->find_gamma_channel( ux ) );
    const int nbin = (upper_channel > lower_channel) ? (upper_channel - lower_channel) : 1;

    // Judge the fit by the statistic it minimized; a likelihood fit generally has a slightly worse
    //  chi2 than the chi2 fit it may be refining.
    const double pre_stat = roi_fit_statistic( inpeaks, data, lower_channel, upper_channel, by_likelihood ) / nbin;
    const double post_stat = roi_fit_statistic( answer, data, lower_channel, upper_channel, by_likelihood ) / nbin;
    const double tolerance = by_likelihood ? (sm_refit_deviance_tolerance / nbin) : 0.0;


    if( (pre_stat + tolerance) < post_stat )
    {
      answer.clear();

      if( true )
      {
        const double ncausalitysigma = 0.0;
        const double stat_threshold  = 0.0;
        const double hypothesis_threshold = 0.0;

        // The caller's choice of fit statistic and uncertainties holds for this refit too.
        Wt::WFlags<PeakFitLM::PeakFitLMOptions> refit_options
                                        = PeakFitLM::PeakFitLMOptions::MediumRefinementOnly;
        for( const PeakFitLMOptions opt : { PeakFitLMOptions::ConditionalAreaUncertainties,
                                            PeakFitLMOptions::NoSparseDataLikelihood,
                                            PeakFitLMOptions::ForcePoissonLikelihood } )
        {
          if( fit_options.test( opt ) )
            refit_options |= opt;
        }
        const vector<shared_ptr<const PeakDef>> refit_peaks
                     = fit_peaks_in_range_LM( lx, ux, ncausalitysigma, stat_threshold, hypothesis_threshold,
                                             inpeaks, data, refit_options, det_type );


        if( refit_peaks.size() == inpeaks.size() )
        {
          const double refit_stat = roi_fit_statistic( refit_peaks, data, lower_channel,
                                                       upper_channel, by_likelihood ) / nbin;
          cout << "refitPeaksThatShareROI_LM: refit_stat=" << refit_stat << ", pre_stat=" << pre_stat << ", post_stat=" << post_stat << endl;
          if( refit_stat <= (pre_stat + tolerance) )
          {
            //cout << "Using re-fit peaks!" << endl;
            answer = refit_peaks;
          }
        }//if( output_peak.size() == inpeaks.size() )
      }

      if( answer.empty() )
        return answer;
    }//if( (pre_stat + tolerance) < post_stat )


    for( const shared_ptr<const PeakDef> &peak : answer )
    {
      const double mean = peak->mean();
      const double sigma = peak->sigma();
      const double data_counts = SpecUtils::gamma_integral( data, mean-0.5*sigma, mean+0.5*sigma );
      const double area = peak->gauss_integral( mean-0.5*sigma, mean+0.5*sigma );

      if( area < 5.0 || area < sqrt(data_counts) )
      {
        answer.clear();
        break;
      }//if( area < 5.0 || area < sqrt(data) )

      if( answer.size() == 1 )
      {
        static int ntimes = 0;
        if( ntimes++ < 3 )
          cerr << "refitPeaksThatShareROI_LM: Should perform some additional chi2 checks for the "
          "case answer.first.size() == 1" << endl;
      }//if( answer.first.size() == 1 )
    }//for( const PeakDefShrdPtr peak : answer.first )

  }catch( std::exception & )
  {
    answer.clear();
    cerr << "refitPeaksThatShareROI_LM: failed to find a better fit" << endl;
  }

  return answer;
}//refitPeaksThatShareROI_LM(...)


/** True when consecutive reweighted (IRLS) passes agree closely enough to stop: every peak's area
 moved by less than 1% of its uncertainty, and its mean and width by less than 0.1% of its width.
 */
static bool irls_pass_converged( const vector<PeakDef> &prev, const vector<PeakDef> &current )
{
  if( prev.size() != current.size() )
    return false;

  for( size_t i = 0; i < current.size(); ++i )
  {
    const PeakDef &a = prev[i], &b = current[i];
    const double sigma = std::max( b.sigma(), 1.0E-9 );
    const double area_tol = 0.01 * std::max( b.amplitudeUncert(), 1.0E-6*fabs(b.amplitude()) );
    if( !(fabs(a.amplitude() - b.amplitude()) <= area_tol)
       || !(fabs(a.mean() - b.mean()) <= 1.0E-3*sigma)
       || !(fabs(a.sigma() - b.sigma()) <= 1.0E-3*sigma) )
      return false;
  }

  return true;
}//irls_pass_converged(...)


/** Internal helper: runs Ceres Levenberg-Marquardt solve on a PeakFitDiffCostFunction.

 For the chi2 objective, ROIs to be fit by Poisson likelihood (the sparse ones by default, every one
 with `ForcePoissonLikelihood`) are then refit by IRLS: the problem is re-solved, warm-started and
 with the same parameter bounds, with each channel's variance frozen at the previous solution's
 model, until the peaks stop changing.  The refit only refines the chi2 fit: if it fails, the chi2
 solution is kept.  (In a joint fit of several ROIs with shared skew, refitting the sparse ROIs also
 moves the shared parameters, so the other ROIs change a little too.)

 Returns the fitted peaks, and the raw fit parameter/uncertainty/covariance data.
 Throws on failure.  If `cancel_flag` is set during the solve, returns with `cancelled` set.
 */
static CeresFitResult run_ceres_fit( PeakFitDiffCostFunction &cost_functor, const size_t total_num_peaks,
                                     const std::shared_ptr<const std::atomic<bool>> &cancel_flag )
{
  const size_t num_fit_pars = cost_functor.number_parameters();
  const PeakFitDiffCostFunction::ProblemSetup prob_setup = cost_functor.get_problem_setup();

  // Which ROIs get reweighted (IRLS) after the chi2 solve: all of them with `ForcePoissonLikelihood`;
  //  by default, the sparse ones (decided from the chi2 solution).  None for the all-in-Ceres
  //  objective, which already is the likelihood fit.
  const bool chi2_objective = (cost_functor.m_objective == FitObjective::NeymanChi2);
  const bool forced_irls = chi2_objective
                           && cost_functor.m_options.test( PeakFitLMOptions::ForcePoissonLikelihood );
  const bool sparse_auto = chi2_objective && !forced_irls
                           && !cost_functor.m_options.test( PeakFitLMOptions::NoSparseDataLikelihood );
  const auto rois_to_reweight = [&]( const double * const pars ) -> vector<bool> {
    if( forced_irls )
      return vector<bool>( cost_functor.m_rois.size(), true );
    if( sparse_auto )
      return cost_functor.sparse_rois( pars );
    return vector<bool>( cost_functor.m_rois.size(), false );
  };

  assert( prob_setup.m_parameters.size() == num_fit_pars );
  assert( prob_setup.m_lower_bounds.size() == num_fit_pars );
  assert( prob_setup.m_upper_bounds.size() == num_fit_pars );

  if( prob_setup.m_constant_parameters.size() >= num_fit_pars )
  {
    // Only the linear solve varies anything, so the reweighting is just iterating it.  A pass that
    //  raises `irls_merit(...)` is undone, and ends it; a failure keeps the chi2 solution.
    const double * const pars = prob_setup.m_parameters.data();
    bool by_likelihood = false;
    try
    {
      const vector<bool> reweight = rois_to_reweight( pars );
      if( sparse_auto )
        tl_fit_objective_diagnostics.sparse_rois += std::count( begin(reweight), end(reweight), true );
      if( std::find( begin(reweight), end(reweight), true ) != end(reweight) )
      {
        vector<vector<double>> models = cost_functor.roi_models( pars );  //of the current solution
        vector<vector<double>> solution_variances;  //the current solution was solved with; empty for chi2
        double merit = cost_functor.irls_merit( models, reweight );
        bool undone = false;
        for( size_t pass = 1; pass <= sm_max_irls_passes; ++pass )
        {
          cost_functor.set_irls_variances( models, reweight );
          vector<vector<double>> new_models = cost_functor.roi_models( pars );
          const double new_merit = cost_functor.irls_merit( new_models, reweight );
          if( new_merit > (merit + sm_irls_deviance_tolerance) )
          {
            if( solution_variances.empty() )
              cost_functor.clear_irls_variances();
            else
              cost_functor.set_irls_variances( solution_variances, reweight );
            undone = true;
            tl_fit_objective_diagnostics.irls_converged = false;
            break;
          }

          double max_rel_change = 0.0;
          for( size_t r = 0; r < models.size(); ++r )
            for( size_t i = 0; i < models[r].size(); ++i )
              max_rel_change = std::max( max_rel_change, fabs(new_models[r][i] - models[r][i])
                                                         / std::max( 1.0, fabs(models[r][i]) ) );
          solution_variances = std::move( models );
          models = std::move( new_models );
          merit = new_merit;
          tl_fit_objective_diagnostics.irls_passes += 1;
          tl_fit_objective_diagnostics.irls_converged = (max_rel_change < 1.0E-6);
          if( tl_fit_objective_diagnostics.irls_converged )
            break;
        }//for( size_t pass = 1; pass <= sm_max_irls_passes; ++pass )

        if( !undone )
          cost_functor.set_irls_variances( models, reweight );
        by_likelihood = !undone || !solution_variances.empty();
      }//if( any ROI to reweight )
    }catch( std::exception & )
    {
      cost_functor.clear_irls_variances();
      by_likelihood = false;
      tl_fit_objective_diagnostics.fell_back_to_chi2 = true;
    }//try / catch

    CeresFitResult result = evaluate_without_free_parameters( cost_functor, prob_setup );
    result.by_likelihood = by_likelihood;
    return result;
  }//if( no free parameters )

  auto cost_function = new ceres::DynamicAutoDiffCostFunction<PeakFitDiffCostFunction,8>(
    &cost_functor, ceres::Ownership::DO_NOT_TAKE_OWNERSHIP );

  cost_function->AddParameterBlock( static_cast<int>(num_fit_pars) );
  cost_function->SetNumResiduals( static_cast<int>(cost_functor.number_residuals()) );

  vector<double> parameters = prob_setup.m_parameters;
  double * const pars = parameters.data();

  std::unique_ptr<ceres::Problem> problem = make_unique<ceres::Problem>();

  // A brief look at a dataset of ~4k HPGe spectra with known truth-value peaks areas shows
  //  that using a loss function doesnt seem improve outcomes, when measured by sucessful
  //  peak fits, and by comparison of fit to truth peak areas.
  ceres::LossFunction *lossfcn = nullptr;
  problem->AddResidualBlock( cost_function, lossfcn, pars ); //Note: problem takes ownership of `cost_function`

  if( !prob_setup.m_constant_parameters.empty() )
  {
    ceres::Manifold *subset_manifold = new ceres::SubsetManifold(
      static_cast<int>(num_fit_pars), prob_setup.m_constant_parameters );
    problem->SetManifold( pars, subset_manifold );
  }

  for( size_t i = 0; i < num_fit_pars; ++i )
  {
    if( prob_setup.m_lower_bounds[i].has_value() )
    {
      assert( *prob_setup.m_lower_bounds[i] <= parameters[i] );
      problem->SetParameterLowerBound( pars, static_cast<int>(i), *prob_setup.m_lower_bounds[i] );
    }

    if( prob_setup.m_upper_bounds[i].has_value() )
    {
      assert( *prob_setup.m_upper_bounds[i] >= parameters[i] );
      problem->SetParameterUpperBound( pars, static_cast<int>(i), *prob_setup.m_upper_bounds[i] );
    }
  }//for( size_t i = 0; i < num_fit_pars; ++i )


  ceres::Solver::Options options;
  options.linear_solver_type = ceres::DENSE_QR;
  options.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH
  options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG
  options.use_nonmonotonic_steps = true;
  options.max_consecutive_nonmonotonic_steps = 10;

  // Initial trust region was very coursely optimized using a fit over ~4k HPGe spectra, using the median peak area
  //  error from truth value (whose value was 1.65*sqrt(area)); but also the average error was also about optimal
  //  (value 3.12*sqrt(area)).
  options.initial_trust_region_radius = 35 * (num_fit_pars - prob_setup.m_constant_parameters.size());
  options.max_trust_region_radius = 1e16;

  // Minimizer terminates when the trust region radius becomes smaller than this value.
  options.min_trust_region_radius = 1e-32;
  // Lower bound for the relative decrease before a step is accepted.
  options.min_relative_decrease = 1e-3; //pseudo optimized based on success rate of fitting peaks - but unknown effect on accuracy of fits.

  // It looks like it usually takes well below 100 calls to solve most problems, but a tail up to ~500, and then rare
  //  (pathelogical?) cases that can take >50k
  //  However, for these pathological cases, if we just try again starting from where it left off, it will solve
  //  it super quick - I guess because there will be more "momentum" to find the minimum, or something, not totally
  //  sure.
  //  We could probably reduce this number of iterations to ~100 - but the effects or impact of this hasnt been
  //  evaluated.
  //  Also, have not checked dependence on number of peaks at all
  const int max_num_iterations = (total_num_peaks > 2) ? 1000 : 500;
  options.max_num_iterations = max_num_iterations;

  // Wall-clock cap is a pathological-case backstop only; the deterministic terminator is
  //  max_num_iterations.  A cap that trips during a normal fit stops after a load-dependent
  //  iteration count => run-to-run nondeterminism, so keep it large.  [determinism fix 2026-07-19]
  options.max_solver_time_in_seconds = 1200.0;
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
  options.minimizer_progress_to_stdout = true;
  options.logging_type = ceres::PER_MINIMIZER_ITERATION;
#else
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;
#endif

  /** Termination when `(new_cost - old_cost) < function_tolerance * old_cost`
   Default 1e-9;  Setting this to 1e-6 from 1e-9 is a ~80x speedup, and slgiht success fitting more peaks.
   Values from 1E-5 to 1E-9 dont seem to have a big effect on accuracy, but 1E-7 had best results for a HPGe dataset

   Value      AvrgAccuracyNSigma   MedianAccuracyNSigma   CPU-Time
   1e-5      3.25068                           1.66096                             465.296 s
   1e-6      3.20428                           1.69388                             581.317 s
   1e-7      3.16456                           1.6934                               1864.5 s
   1e-8      3.29325                           1.6936                               16833.9 s
   */
  options.function_tolerance = 1e-7;
  options.parameter_tolerance = 1e-11; //Default value is 1e-8.  Using 1e-11, so its usually the function tolerance that terminates things.
  options.num_threads = 1; //Probably wont have much/any effect

  // Default value of `max_num_consecutive_invalid_steps` is 5, however, if we are re-fitting a peak who already has
  //  near perfect values, occasionally the fit fails with the devault values - I guess the trust-region is just so
  //  far off (but I dont really know - this is just a guess), but if we increase this value to 10, the fit seems
  //  to be sucessful, for a limited number of test cases.
  options.max_num_consecutive_invalid_steps = 10;

  // Cooperative cancellation: if a `cancel_flag` was provided (e.g., by the automated peak search
  //  during application shutdown), an iteration callback makes the solver bail within one LM step
  //  rather than running out the full `max_solver_time_in_seconds`.  Ceres reports `USER_FAILURE`.
  std::unique_ptr<CancelIterationCallback> cancel_callback;
  if( cancel_flag )
  {
    cancel_callback.reset( new CancelIterationCallback( cancel_flag ) );
    options.callbacks.push_back( cancel_callback.get() );
  }

  const auto was_cancelled = [&cancel_flag]() -> bool {
    return cancel_flag && cancel_flag->load( std::memory_order_relaxed );
  };

  // Solves from the current `parameters`, retrying from where it stopped if it did not converge, and
  //  then with numerical differentiation if that failed too (which replaces `problem`).  Returns
  //  false if cancelled; throws if the solve fails.
  const auto solve = [&]() -> bool {
    options.max_num_iterations = max_num_iterations;

    ceres::Solver::Summary summary;
    ceres::Solve( options, problem.get(), &summary );

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    std::cout << summary.FullReport() << "\n";
    cout << "run_ceres_fit: Took " << cost_functor.m_ncalls.load() << " calls to solve." << endl;
#endif

    // Cancellation (`SOLVER_ABORT`, reported as `USER_FAILURE`) short-circuits the retries, which
    //  would just abort again.
    if( was_cancelled() )
      return false;

    string failure_reason;

    switch( summary.termination_type )
    {
      case ceres::CONVERGENCE:
      case ceres::USER_SUCCESS:
        break;

      case ceres::NO_CONVERGENCE:
      {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "run_ceres_fit: NO_CONVERGENCE:\n" << summary.FullReport() << endl;
#endif
        // We will give it another go - for most cases it seems the solution will now be succesffuly found really
        //  quickly.  The guess is this is from momentum, or whatever, but not really sure.
        options.max_num_iterations = 5000;
        cost_functor.m_ncalls = 0;
        summary = ceres::Solver::Summary();
        ceres::Solve( options, problem.get(), &summary );

        if( was_cancelled() )
          return false;

        if( (summary.termination_type == ceres::CONVERGENCE)
           || (summary.termination_type == ceres::USER_SUCCESS) )
        {
          break;
        }

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "run_ceres_fit: Retry failed with " << cost_functor.m_ncalls.load() << " additional calls" << endl;
#endif
        failure_reason = "The L-M ceres::Solver solving failed - NO_CONVERGENCE.";
        [[fallthrough]];
      }//case ceres::NO_CONVERGENCE:

      case ceres::FAILURE:
      {
        if( failure_reason.empty() )
          failure_reason = "The L-M ceres::Solver solving failed - FAILURE.";

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "run_ceres_fit: FAILURE:\n" << summary.FullReport() << endl;
#endif
        // We will re-try with numerical differntiation - the one case that the author has observed this being
        //  necassary is with peak-refits, where all peak paramaeters are already about perfect, very rarely, it seems
        //  to throw off trust-region searching, or something - however, numerical differentiating seems to work for
        //  these (rare!) cases that this is observed.
        auto numeric_cost_fnct = new ceres::DynamicNumericDiffCostFunction<PeakFitDiffCostFunction>(
          &cost_functor, ceres::Ownership::DO_NOT_TAKE_OWNERSHIP );
        numeric_cost_fnct->AddParameterBlock( static_cast<int>(num_fit_pars) );
        numeric_cost_fnct->SetNumResiduals( static_cast<int>(cost_functor.number_residuals()) );

        std::unique_ptr<ceres::Problem> numerical_problem = make_unique<ceres::Problem>();
        numerical_problem->AddResidualBlock( numeric_cost_fnct, lossfcn, pars );

        if( !prob_setup.m_constant_parameters.empty() )
        {
          ceres::Manifold *subset_manifold = new ceres::SubsetManifold(
            static_cast<int>(num_fit_pars), prob_setup.m_constant_parameters );
          numerical_problem->SetManifold( pars, subset_manifold );
        }

        for( size_t i = 0; i < num_fit_pars; ++i )
        {
          if( prob_setup.m_lower_bounds[i].has_value() )
            numerical_problem->SetParameterLowerBound( pars, static_cast<int>(i), *prob_setup.m_lower_bounds[i] );
          if( prob_setup.m_upper_bounds[i].has_value() )
            numerical_problem->SetParameterUpperBound( pars, static_cast<int>(i), *prob_setup.m_upper_bounds[i] );
        }//for( size_t i = 0; i < num_fit_pars; ++i )

        cost_functor.m_ncalls = 0;
        summary = ceres::Solver::Summary();
        ceres::Solve( options, numerical_problem.get(), &summary );

        if( was_cancelled() )
          return false;

        switch( summary.termination_type )
        {
          case ceres::CONVERGENCE:
          case ceres::USER_SUCCESS:
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
            cout << "run_ceres_fit: Solved using numerical differentiation." << endl;
#endif
            problem = std::move( numerical_problem );
            break;

          case ceres::NO_CONVERGENCE:
          case ceres::FAILURE:
          case ceres::USER_FAILURE:
          {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
            cerr << "run_ceres_fit: Failed with numerical differentiation too - giving up." << endl;
#endif
            throw runtime_error( failure_reason );
          }
        }//switch( summary.termination_type )

        break;
      }//case ceres::FAILURE:

      case ceres::USER_FAILURE:
      {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "run_ceres_fit: USER_FAILURE:\n" << summary.FullReport() << endl;
#endif
        throw runtime_error( "The L-M ceres::Solver solving failed - USER_FAILURE." );
      }
    }//switch( summary.termination_type )

    return true;
  };//solve lambda

  CeresFitResult result;
  result.num_fit_pars = num_fit_pars;
  result.by_likelihood = !chi2_objective;

  if( !solve() )
  {
    result.cancelled = true;
    return result;
  }

  // Refit the ROIs that call for it by IRLS.  This only refines the chi2 fit, so if anything in it
  //  fails, the chi2 solution is kept.
  const vector<double> chi2_parameters = parameters;
  try
  {
    const vector<bool> reweight = rois_to_reweight( pars );
    if( sparse_auto )
      tl_fit_objective_diagnostics.sparse_rois += std::count( begin(reweight), end(reweight), true );
    if( std::find( begin(reweight), end(reweight), true ) != end(reweight) )
    {
      // Reweight: freeze each channel's variance at the current model and solve again, until the
      //  peaks stop changing.  The parameter bounds stay those of the input peaks.  ROIs not being
      //  reweighted keep their chi2 weighting throughout.
      //
      //  The reweighted objective has the same gradient as `irls_merit(...)` (for a single ROI, the
      //  Poisson deviance) at the point its variances were frozen at, so each pass's step is a
      //  descent direction for it - but a full step can overshoot (with near-empty channels it can
      //  wander off to a far worse solution).  So we backtrack along the step until the merit
      //  drops, and stop if it will not.
      vector<vector<double>> models;
      vector<PeakDef> prev_peaks = cost_functor.parametersToPeaks<PeakDef,double>( pars, nullptr, nullptr,
                                                                                  nullptr, 0, &models );
      double best_deviance = cost_functor.irls_merit( models, reweight );

      for( size_t pass = 1; pass <= sm_max_irls_passes; ++pass )
      {
        const vector<double> start_parameters = parameters;
        cost_functor.set_irls_variances( models, reweight );
        cost_functor.m_ncalls = 0;

        bool solved = false;
        try
        {
          solved = solve();
          if( !solved )
          {
            result.cancelled = true;
            return result;
          }
        }catch( std::exception & )
        {
          // Keep the best pass so far.
          std::copy( begin(start_parameters), end(start_parameters), begin(parameters) );
          tl_fit_objective_diagnostics.irls_converged = false;
          break;
        }

        tl_fit_objective_diagnostics.irls_passes += 1;

        // Ceres stops each solve at a relative cost change of `function_tolerance`, so the deviance of
        //  an already-converged point can come back a hair higher; allow for that much.
        const double prev_deviance = best_deviance;
        const double accept_tol = options.function_tolerance * std::max( 1.0, fabs(best_deviance) );
        const vector<double> full_step = parameters;
        double full_step_deviance = std::numeric_limits<double>::infinity();
        bool accepted = false;
        for( double frac = 1.0; frac > 0.06; frac *= 0.5 )
        {
          for( size_t i = 0; i < num_fit_pars; ++i )
            parameters[i] = start_parameters[i] + frac*(full_step[i] - start_parameters[i]);

          vector<vector<double>> trial_models = cost_functor.roi_models( pars );
          const double deviance = cost_functor.irls_merit( trial_models, reweight );
          if( frac == 1.0 )
            full_step_deviance = deviance;

          if( deviance <= (best_deviance + accept_tol) )
          {
            accepted = true;
            best_deviance = std::min( deviance, best_deviance );
            models = std::move( trial_models );
            break;
          }
        }//for( backtracking along the step )

        if( !accepted )
        {
          // Back to where this pass started.  If even the full step barely changed the deviance, we
          //  were already at the minimum.
          std::copy( begin(start_parameters), end(start_parameters), begin(parameters) );
          tl_fit_objective_diagnostics.irls_converged
                          = ((full_step_deviance - prev_deviance) < sm_irls_deviance_tolerance);
          break;
        }

        vector<PeakDef> peaks = cost_functor.parametersToPeaks<PeakDef,double>( pars, nullptr, nullptr,
                                                                               nullptr, 0, nullptr );
        const bool converged = irls_pass_converged( prev_peaks, peaks )
                               || ((prev_deviance - best_deviance) < sm_irls_deviance_tolerance);
        prev_peaks = std::move( peaks );

        tl_fit_objective_diagnostics.irls_converged = converged;
        if( converged )
          break;
      }//for( size_t pass = 1; pass <= sm_max_irls_passes; ++pass )

      // Freeze the variances at the final model, so the covariance below is the expected Fisher
      //  information at the reported point (and the linear solve takes its last reweighting step).
      cost_functor.set_irls_variances( cost_functor.roi_models( pars ), reweight );
      result.by_likelihood = true;
    }//if( any ROI to reweight )
  }catch( std::exception & )
  {
    std::copy( begin(chi2_parameters), end(chi2_parameters), begin(parameters) );
    cost_functor.clear_irls_variances();
    result.by_likelihood = false;
    tl_fit_objective_diagnostics.fell_back_to_chi2 = true;
    tl_fit_objective_diagnostics.irls_converged = false;
  }//try / catch( IRLS refit )


  // Compute covariance
  cost_functor.m_ncalls = 0;

  ceres::Covariance::Options cov_options;
  cov_options.algorithm_type = ceres::CovarianceAlgorithmType::DENSE_SVD; //SPARSE_QR;

  // Some terms we are fitting may not matter much, so when computing the inverse of J'J, we will opt
  //  to drop all terms where the eigen value of that term, divided by the maximum eigenvalue is
  //  less than the condition number.
  //  TODO: Currently leaving condition number at default 1e-14 - but what value is reasonable should be investigated
  cov_options.null_space_rank = -1;
  cov_options.min_reciprocal_condition_number = 1e-14;

  vector<double> uncertainties( num_fit_pars, 0.0 );
  double *uncertainties_ptr = uncertainties.data();
  vector<double> row_major_covariance;

  ceres::Covariance covariance( cov_options );
  vector<pair<const double*, const double*>> covariance_blocks;
  covariance_blocks.push_back( make_pair( pars, pars ) );

  // For the deviance residuals, the covariance comes from Fisher residuals at the fit's model; see
  //  `PeakFitDiffCostFunction::m_fisher_variances`.
  const bool fisher_covariance = (cost_functor.m_objective == FitObjective::PoissonAllInCeres);
  if( fisher_covariance )
    cost_functor.set_fisher_variances( cost_functor.roi_models( pars ) );

  const bool covariance_ok = covariance.Compute( covariance_blocks, problem.get() );

  if( fisher_covariance )
    cost_functor.set_fisher_variances( {} );

  if( !covariance_ok )
  {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
    cerr << "run_ceres_fit: Failed to compute covariance!" << endl;
#endif
    uncertainties_ptr = nullptr;
  }
  else
  {
    row_major_covariance.resize( num_fit_pars * num_fit_pars );
    const vector<const double *> const_par_blocks( 1, pars );

    const bool success = covariance.GetCovarianceMatrix( const_par_blocks, row_major_covariance.data() );
    assert( success );
    if( success )
    {
      for( size_t i = 0; i < num_fit_pars; ++i )
      {
        if( row_major_covariance[i * num_fit_pars + i] > 0.0 )
          uncertainties[i] = sqrt( row_major_covariance[i * num_fit_pars + i] );
      }
    }
    else
    {
      uncertainties_ptr = nullptr;
      row_major_covariance.clear();
    }
  }//if( failed covariance ) / else

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
  cout << "run_ceres_fit: Took " << cost_functor.m_ncalls.load() << " calls to get covariance." << endl;
#endif

  const double *cov_ptr = row_major_covariance.empty() ? nullptr : row_major_covariance.data();
  vector<double> residuals( cost_functor.number_residuals(), 0.0 );
  result.final_peaks = cost_functor.parametersToPeaks<PeakDef,double>(
    parameters.data(), uncertainties_ptr, residuals.data(), cov_ptr, num_fit_pars );

  if( cov_ptr && !cost_functor.m_options.test( PeakFitLMOptions::ConditionalAreaUncertainties ) )
  {
    try
    {
      cost_functor.add_nonlinear_uncertainty_to_linear_pars( parameters.data(), row_major_covariance,
                                                             prob_setup, result.final_peaks );
    }catch( std::exception & )
    {
      // Keep the conditional uncertainties; the peaks are only changed once nothing can throw.
    }
  }

  result.parameters = std::move( parameters );
  result.uncertainties = std::move( uncertainties );
  result.row_major_covariance = std::move( row_major_covariance );

  return result;
}//run_ceres_fit(...)


FitPeaksResults fit_peaks_in_spectrum_LM( const vector<shared_ptr<const PeakDef>> input_peaks,
                                          shared_ptr<const SpecUtils::Measurement> data,
                                          const double stat_threshold,
                                          const double hypothesis_threshold,
                                          const std::optional<PeakFitUtils::CoarseResolutionType> resolution_type,
                                          const std::optional<PeakDef::SkewType> skew_type,
                                          const Wt::WFlags<PeakFitLMOptions> fit_options,
                                          const std::function<bool(const PeakDef &,const PeakDef &)>
                                            &may_remove_close_pair ) throw()
{
  FitPeaksResults results;

  try
  {
    if( !data || !data->gamma_counts() || data->gamma_counts()->empty() )
    {
      results.status = FitPeaksResults::FitPeaksResultsStatus::Failure;
      results.error_message = "fit_peaks_in_spectrum_LM: invalid spectrum data.";
      return results;
    }

    if( input_peaks.empty() )
    {
      results.status = FitPeaksResults::FitPeaksResultsStatus::Failure;
      results.error_message = "fit_peaks_in_spectrum_LM: no input peaks.";
      return results;
    }

    // Separate data-defined (non-Gaussian) peaks from Gaussian peaks
    vector<shared_ptr<const PeakDef>> gauss_peaks;
    vector<shared_ptr<const PeakDef>> datadefined_peaks;
    gauss_peaks.reserve( input_peaks.size() );

    for( const shared_ptr<const PeakDef> &p : input_peaks )
    {
      if( p->gausPeak() )
        gauss_peaks.push_back( p );
      else
        datadefined_peaks.push_back( p );
    }

    if( gauss_peaks.empty() )
    {
      // No Gaussian peaks to fit - just return data-defined peaks as-is
      results.status = FitPeaksResults::FitPeaksResultsStatus::Success;
      results.fit_peaks = datadefined_peaks;
      return results;
    }

    // Determine detector resolution type
    const PeakFitUtils::CoarseResolutionType res_type = resolution_type.has_value()
      ? *resolution_type
      : PeakFitUtils::coarse_resolution_from_peaks( gauss_peaks );

    // Deep copy all Gaussian peaks so we dont modify the inputs
    local_unique_copy_continuum( gauss_peaks );

    vector<shared_ptr<const PeakDef>> all_fit_peaks;

    if( !skew_type.has_value() )
    {
      // Mode A: skew_type is nullopt - fit each ROI independently
      //  Group peaks by shared continuum
      map<shared_ptr<const PeakContinuum>, vector<shared_ptr<const PeakDef>>> cont_to_peaks;
      for( const shared_ptr<const PeakDef> &p : gauss_peaks )
        cont_to_peaks[p->continuum()].push_back( p );

      for( const auto &kv : cont_to_peaks )
      {
        try
        {
          const vector<shared_ptr<const PeakDef>> roi_result
            = fit_peaks_in_roi_LM( kv.second, data, res_type, fit_options );
          all_fit_peaks.insert( end(all_fit_peaks), begin(roi_result), end(roi_result) );
        }
        catch( std::exception &e )
        {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
          cerr << "fit_peaks_in_spectrum_LM: ROI fit failed: " << e.what() << endl;
#endif
          // If a ROI fails, keep the original peaks for that ROI
          all_fit_peaks.insert( end(all_fit_peaks), begin(kv.second), end(kv.second) );
        }
      }//for( each ROI group )
    }
    else
    {
      // Mode B: skew_type is specified - fit all ROIs simultaneously with shared/related skew
      const PeakDef::SkewType target_skew = *skew_type;

      // Set skew type on all peaks (with VoigtPlusBortel exception)
      vector<shared_ptr<const PeakDef>> peaks_to_fit;
      peaks_to_fit.reserve( gauss_peaks.size() );

      for( const shared_ptr<const PeakDef> &p : gauss_peaks )
      {
        shared_ptr<PeakDef> newpeak = make_shared<PeakDef>( *p );

        const bool keep_voigt = (newpeak->skewType() == PeakDef::SkewType::VoigtPlusBortel)
          && ((target_skew == PeakDef::SkewType::Bortel) || (target_skew == PeakDef::SkewType::GaussPlusBortel));

        if( !keep_voigt && (newpeak->skewType() != target_skew) )
        {
          newpeak->setSkewType( target_skew );

          // Set default starting values and enable fitting for skew parameters
          const size_t num_skew_pars = PeakDef::num_skew_parameters( target_skew );
          for( size_t i = 0; i < num_skew_pars; ++i )
          {
            const PeakDef::CoefficientType ct
              = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );
            double lower, upper, start, dx;
            PeakDef::skew_parameter_range( target_skew, ct, lower, upper, start, dx );
            newpeak->set_coefficient( start, ct );
            newpeak->setFitFor( ct, true );
          }
        }
        // else if peak already had the target skew type, preserve its fitFor settings

        peaks_to_fit.push_back( newpeak );
      }//for( each gauss peak )

      try
      {
        // PeakFitDiffCostFunction handles VoigtPlusBortel exception and multi-ROI grouping internally
        PeakFitDiffCostFunction cost_functor( data, peaks_to_fit, 0, 0, 0,
                                              target_skew, res_type, fit_options );

        CeresFitResult ceres_result = run_ceres_fit( cost_functor, peaks_to_fit.size() );

        // Convert to shared_ptr
        all_fit_peaks.reserve( ceres_result.final_peaks.size() );
        for( PeakDef &peak : ceres_result.final_peaks )
          all_fit_peaks.push_back( make_shared<PeakDef>( std::move(peak) ) );

        // Extract SkewRelation if applicable
        const size_t num_skew_pars = PeakDef::num_skew_parameters( target_skew );
        const bool should_have_skew_relation = (cost_functor.m_rois.size() > 1)
          && !fit_options.test( PeakFitLMOptions::IndependentSkewValues )
          && (num_skew_pars > 0);

        if( should_have_skew_relation )
        {
          FitPeaksResults::SkewRelation skew_rel;
          skew_rel.skew_type = target_skew;
          skew_rel.energy_range = std::make_pair( cost_functor.m_skew_anchor_lower_energy,
                                                   cost_functor.m_skew_anchor_upper_energy );

          size_t upper_offset = num_skew_pars;
          for( size_t i = 0; i < num_skew_pars; ++i )
          {
            const PeakDef::CoefficientType ct
              = PeakDef::CoefficientType( static_cast<int>(PeakDef::SkewPar0) + static_cast<int>(i) );

            if( PeakDef::is_energy_dependent( target_skew, ct ) && cost_functor.m_fit_skew_energy_dependence )
            {
              skew_rel.energy_dependent_skew_pars[i]
                = std::make_pair( ceres_result.parameters[i], ceres_result.parameters[upper_offset] );
              upper_offset += 1;
            }
            else
            {
              skew_rel.non_energy_dependent_skew_pars[i]
                = std::make_pair( ceres_result.parameters[i], ceres_result.uncertainties[i] );
            }
          }//for( each skew parameter )

          results.skew_relation = skew_rel;
        }//if( should_have_skew_relation )
      }
      catch( std::exception &e )
      {
#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO )
        cerr << "fit_peaks_in_spectrum_LM: multi-ROI fit failed: " << e.what() << endl;
#endif
        results.status = FitPeaksResults::FitPeaksResultsStatus::Failure;
        results.error_message = string( "fit_peaks_in_spectrum_LM: " ) + e.what();
        return results;
      }
    }//if( !skew_type.has_value() ) / else


    // Apply significance testing - mirroring fit_peaks_LM
    for( size_t i = 0; i < all_fit_peaks.size(); )
    {
      const shared_ptr<const PeakDef> &peak = all_fit_peaks[i];

      // Dont remove peaks whose amplitudes we arent fitting
      if( !peak->fitFor( PeakDef::GaussAmplitude ) )
      {
        ++i;
        continue;
      }

      const double num_sigma = (peak->amplitudeUncert() > 0.0)
        ? (peak->amplitude() / peak->amplitudeUncert()) : 999.0;
      const bool is_sig = (stat_threshold <= 0.0) || (num_sigma >= stat_threshold);

      const double dummy_stat_thresh = 0.0;
      const bool significant = chi2_significance_test( *peak, dummy_stat_thresh, hypothesis_threshold, {}, data );

      if( !is_sig || !significant )
      {
        results.lost_peaks.push_back( peak );
        all_fit_peaks.erase( all_fit_peaks.begin() + i );
      }
      else
      {
        ++i;
      }
    }//for( significance testing )

    // Remove peaks whose means are within 1 sigma of each other
    std::sort( begin(all_fit_peaks), end(all_fit_peaks), &PeakDef::lessThanByMeanShrdPtr );

    for( size_t i = 1; i < all_fit_peaks.size(); ++i )
    {
      const shared_ptr<const PeakDef> &this_peak = all_fit_peaks[i - 1];
      const shared_ptr<const PeakDef> &next_peak = all_fit_peaks[i];

      if( !this_peak->fitFor( PeakDef::GaussAmplitude ) )
        continue;

      const double min_sigma = std::min( this_peak->sigma(), next_peak->sigma() );
      const double mean_diff = next_peak->mean() - this_peak->mean();

      if( ((mean_diff / min_sigma) < 1.0)
          && (!may_remove_close_pair || may_remove_close_pair(*this_peak, *next_peak)) )
      {
        // Remove the peak with the worse chi2
        if( this_peak->chi2dof() > next_peak->chi2dof() )
        {
          results.lost_peaks.push_back( this_peak );
          all_fit_peaks.erase( all_fit_peaks.begin() + static_cast<ptrdiff_t>(i - 1) );
        }
        else
        {
          results.lost_peaks.push_back( next_peak );
          all_fit_peaks.erase( all_fit_peaks.begin() + static_cast<ptrdiff_t>(i) );
        }
        i = (i > 1) ? (i - 1) : 0;
      }
    }//for( removing duplicate peaks )

    // Add back data-defined peaks
    if( !datadefined_peaks.empty() )
    {
      all_fit_peaks.insert( end(all_fit_peaks), begin(datadefined_peaks), end(datadefined_peaks) );
      std::sort( begin(all_fit_peaks), end(all_fit_peaks), &PeakDef::lessThanByMeanShrdPtr );
    }

    results.status = FitPeaksResults::FitPeaksResultsStatus::Success;
    results.fit_peaks = std::move( all_fit_peaks );
  }
  catch( std::exception &e )
  {
    results.status = FitPeaksResults::FitPeaksResultsStatus::Failure;
    results.error_message = string( "fit_peaks_in_spectrum_LM: " ) + e.what();

#if( PRINT_VERBOSE_PEAK_FIT_LM_INFO || !defined(NDEBUG) )
    cerr << results.error_message << endl;
#endif
  }

  return results;
}//fit_peaks_in_spectrum_LM(...)


}//namespace PeakFitLM
