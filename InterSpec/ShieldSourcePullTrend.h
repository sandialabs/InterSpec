#ifndef ShieldSourcePullTrend_h
#define ShieldSourcePullTrend_h
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

#include <memory>
#include <string>
#include <vector>
#include <utility>

// Forward declarations
class DetectorPeakResponse;

namespace GammaInteractionCalc
{
  struct PeakDetail;
  struct PeakResultPlotInfo;
}//namespace GammaInteractionCalc

namespace ShieldingSourceFitCalc
{
  struct SourceFitDef;
  struct ShieldingInfo;
}//namespace ShieldingSourceFitCalc


/** Diagnoses likely shielding-model problems from the per-peak pulls ("chi" values) of an
 activity/shielding fit.

 The decision uses profile likelihood-ratio tests along the PHYSICAL directions a wrong shielding
 model can be off in, evaluated from what the fit already reported for each peak (its observed and
 model counts, the per-nuclide model contributions, and the model shielding optical depth):
  - amount: scale the model shieldings optical depth - or, when the model has no shielding, add a
    reference layer (Z 7, 13, 26, 50 or 82) - re-solving every nuclides activity;
  - atomic number: free the effective Z and areal density of a single-shield model;
  - "other": a shared smooth (quadratic in log energy) efficiency multiplier beyond the best of those.
 Each test is the signed square root of the chi2 improvement, divided by a residual scale estimated
 from the fit (never below a prior of 2, since real peak fits scatter 2-5x their statistical
 uncertainty), with a 2% fractional uncertainty floor so a few very strong peaks cannot decide alone.

 A low-order polynomial is still fit to the pulls vs log(energy), but only to draw a descriptive
 trend curve on the chi chart; it no longer decides anything.  (Its uncertainties are inflated by its
 own residual scatter, so a large step-shaped error - the usual signature of wrong shielding - used to
 mask itself; and its curvature pre-empted amount calls.)

 Calibrated by scratch/20261111_pull_trend_study (BASELINE.md): GADRAS inject corpus truth, re-
 attenuated truth and real FitPeaksForNuclides fits, across 19 detectors and 30-1800 s live times.
 */
namespace ShieldSourcePullTrend
{
  /** The transform of energy used as the regressor for the trend fit.
   The default basis is chosen from the calibration study (see implementation notes).
   */
  enum class TrendBasis : int
  {
    Energy,
    LogEnergy,
    InvEnergy,
    /** Optical depth of the model shielding at each energy - i.e., sum of mu*thickness.
     Requires the shielding definitions; falls back to LogEnergy if they are unusable.
     */
    MuModel
  };//enum class TrendBasis


  /** Which trend model ended up statistically supported by the points. */
  enum class TrendModel : int
  {
    None,
    Linear,
    Quadratic
  };//enum class TrendModel


  /** The diagnosis reached from the trend, gated by what the fit was allowed to vary. */
  enum class TrendConclusion : int
  {
    NoConclusion,
    /** The model appears to have less shielding than is actually present. */
    ShieldingAmountTooLow,
    /** The model appears to have more shielding than is actually present. */
    ShieldingAmountTooHigh,
    /** The shielding's effective atomic number appears lower than reality. */
    AtomicNumberTooLow,
    /** The shielding's effective atomic number appears higher than reality. */
    AtomicNumberTooHigh,
    /** A smooth residual energy trend remains that the shielding changes do not explain - the
     distance, geometry, source age, or detector response may be off.
     */
    OtherModelInconsistency,
    /** The model has no shielding, but the peaks clearly want some (see
     TrendResult::suggestedArealDensity and TrendResult::suggestedAtomicNumber).
     */
    ModelLacksShielding
  };//enum class TrendConclusion


  /** How strongly the test statistic supports the conclusion; controls message wording
   ("may be" / "probably" / "strongly suggests").  Tier boundaries come from the
   calibration study.
   */
  enum class TrendConfidence : int
  {
    /** A soft advisory ("possibly"): an atomic-number signal of 2.5-3.5 sigma, only reported when
     the data demonstrably have the power to resolve a Z error.  Only reachable for atomic-number calls.
     */
    Tentative,
    Weak,
    Moderate,
    Strong
  };//enum class TrendConfidence


  /** Summary of what the activity/shielding fit was allowed to vary - used to gate which
   conclusions are logically possible (e.g., if shielding thickness was fit, the fit should
   have found the optimal amount, so only an atomic-number diagnosis makes sense).
   */
  struct FitConfigSummary
  {
    /** Any dimension of a physical (non-generic) shielding was fit (thickness, radius, etc). */
    bool physicalDimensionFit = false;

    /** A generic shielding's atomic number was fit. */
    bool genericAnFit = false;

    /** A generic shielding's areal density was fit. */
    bool genericAdFit = false;

    /** A self-attenuating (intrinsic) or trace source is present; the simple point-source
     attenuation reasoning does not apply, so conclusions are suppressed (curves still drawn).
     */
    bool hasSelfAttenOrTraceSrc = false;

    /** A source age was fit; age changes line ratios in ways the shielding tests do not model, so
     only the "other" conclusion is possible.
     */
    bool anyAgeFit = false;

    /** The source is close to the detector (see #set_summing_info), cascade-summing correction is
     off, and a nuclide has coincident gammas: the residuals may partly be summing, so any message
     carries a caveat.  Filled by #set_summing_info (false otherwise).
     */
    bool possibleCascadeSumming = false;
  };//struct FitConfigSummary


  /** Builds a #FitConfigSummary from the shielding and source definitions of a fit. */
  FitConfigSummary summarize_fit_config(
                      const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                      const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources );

  /** Sets FitConfigSummary::possibleCascadeSumming: true when `distance` (source to detector face)
   is under 15 cm, `summing_corrected` is false, some source nuclide (or its progeny) emits coincident
   gammas, and the detector's maximum summing magnitude is at least ~1% (or cannot be estimated
   because the DRF has no total-efficiency information).  Cheap; never throws.
   */
  void set_summing_info( FitConfigSummary &config, const double distance, const bool summing_corrected,
                         const std::shared_ptr<const DetectorPeakResponse> &drf,
                         const std::vector<ShieldingSourceFitCalc::SourceFitDef> &sources );


  /** The trend curve for one nuclide - the shared shape offset by this nuclides intercept. */
  struct NuclideTrend
  {
    /** Nuclide symbol, as given by `PeakDetail::assignedNuclide`. */
    std::string nuclide;

    /** CSS color of the first peak of this nuclide with a non-empty color; may be empty. */
    std::string cssColor;

    size_t numPeaksUsed = 0;

    /** Energy range, in keV, of this nuclides used peaks; curves are not extrapolated past it. */
    double minEnergy = 0.0, maxEnergy = 0.0;

    double intercept = 0.0, interceptUncert = 0.0;

    /** The fitted trend chi(E), sampled across [minEnergy,maxEnergy]; {energy_keV, chi} pairs.
     Empty if the selected model is TrendModel::None.
     */
    std::vector<std::pair<double,double>> curve;
  };//struct NuclideTrend


  /** Result of the pull-trend diagnosis. */
  struct TrendResult
  {
    TrendBasis basis = TrendBasis::LogEnergy;
    TrendModel model = TrendModel::None;

    /** Shared slope/curvature in the centered, [-1,1]-scaled regressor coordinate, with
     their 1-sigma uncertainties and t-statistics (coefficient/uncertainty).
     Note: because the regressor is normalized, values are comparable across bases.
     */
    double slope = 0.0, slopeUncert = 0.0, slopeT = 0.0;
    double curvature = 0.0, curvatureUncert = 0.0, curvatureT = 0.0;

    /** Sum of squared residuals of the selected trend model, and its degrees of freedom. */
    double chi2 = 0.0;
    unsigned int dof = 0;

    /** Residual variance of the pulls about the trend (chi2/dof).  Under a correct model with
     well-behaved peak fits this is ~1; values >> 1 indicate over-dispersed pulls (2x-Poisson
     peak-fit systematics, or a poorly-converged shielding fit).  The reported t-statistics
     already divide the coefficient uncertainties by sqrt(max(1,redChi2)) so they stay
     calibrated, which both controls the false-positive rate and suppresses garbage from
     poorly-converged fits.
     */
    double redChi2 = 0.0;

    /** Total pull points used in the fit (after nuclide-ambiguity and validity filtering). */
    size_t numPointsUsed = 0;

    /** Profile likelihood-ratio statistics that drive the conclusion (normal-equivalent sigma).
     amountT > 0: the peaks want MORE shielding than the model has (< 0: less).
     atomicNumberT > 0: the peaks want a HIGHER effective atomic number than the model (model Z too
     low); < 0: lower.  otherT >= 0: a smooth residual trend beyond the best physical alternative.
     */
    double amountT = 0.0, atomicNumberT = 0.0, otherT = 0.0;

    /** Whether the likelihood-ratio tests ran (not for self-attenuating or trace sources, or with
     too few peaks); when false the statistics below are 0 and mean "not tested", not "nothing found".
     */
    bool testsRun = false;

    /** Residual scale the statistics were divided by (>= the prior of 2); 0 when the tests did not run. */
    double residualScale = 0.0;

    /** Peaks (after validity filtering) the likelihood-ratio tests used. */
    size_t numPeaksTested = 0;

    /** Whether the model had any shielding with non-zero areal density. */
    bool modelHasShielding = false;

    /** Best-fit total areal density (PhysicalUnits) of the amount alternative - the added layer when
     the model had no shielding, else the scaled model shielding; 0 if not computed.
     */
    double suggestedArealDensity = 0.0;

    /** Best-fit effective atomic number: of the added layer (ModelLacksShielding), or of the free-Z
     alternative (atomic-number conclusions); 0 if not computed.
     */
    double suggestedAtomicNumber = 0.0;

    /** For ModelLacksShielding: whether the reference-Z layers are distinguishable enough to name a
     Z class (otherwise only an iron-equivalent areal density is meaningful).
     */
    bool suggestedAtomicNumberResolved = false;

    /** Copied from FitConfigSummary::possibleCascadeSumming. */
    bool possibleCascadeSumming = false;

    std::vector<NuclideTrend> nuclideTrends;

    /** Atomic-number DISCRIMINATION POWER: the atomicNumberT that a reference effective-Z error (a
     factor 1.5, either direction) WOULD produce for this exact peak set and model shielding, if it
     existed - a design property of the data, not of the observed pulls.  Above ~250 keV attenuation
     is Compton-dominated and nearly Z-independent in shape, so only low-energy, significant peaks
     carry Z information; a small power means a null result is "inconclusive", not "Z consistent".
     0 when not computed (no single shielding in the model, or too few peaks).
     */
    double anDiscriminationT = 0.0;

    /** True when #anDiscriminationT reaches the power gate (2.5): soft (2.5-3.5 sigma) atomic-number
     calls are only made then; firmer ones (>= 3.5 sigma) speak for themselves.
     */
    bool anDiscriminable = false;

    TrendConclusion conclusion = TrendConclusion::NoConclusion;
    TrendConfidence confidence = TrendConfidence::Weak;
  };//struct TrendResult


  /** Diagnoses the fit's per-peak pulls (see the namespace description), and fits the descriptive
   pull-vs-energy curve (per-nuclide intercepts plus shared slope/curvature, in centered/scaled `basis`
   coordinates) that the chi chart draws.

   \param peak_comparisons Per-peak pull info, as from ShieldingSourceChi2Fcn::energy_chi_contributions
   \param peak_details Parallel per-peak details (for nuclide assignment and source fractions);
          must correspond index-for-index with `peak_comparisons` (extra entries in either are ignored)
   \param shieldings The model shieldings (as fit).  Required: the amount and atomic-number tests are
          anchored on them; an empty vector means "the model has no shielding".
   \param config What the fit varied - gates the possible conclusions
   \param basis The regressor transform of energy

   Never throws for degenerate input; returns a result with model==TrendModel::None instead.
   */
  TrendResult fit_pull_trend(
                const std::vector<GammaInteractionCalc::PeakResultPlotInfo> &peak_comparisons,
                const std::vector<GammaInteractionCalc::PeakDetail> &peak_details,
                const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
                const FitConfigSummary &config,
                const TrendBasis basis = TrendBasis::LogEnergy );


  /** Localization key (in ShieldingSourceDisplay.xml) for the conclusion message; the message
   takes the qualifier text (see #confidence_txt_key) as argument {1}, and the hint (see
   #message_hint) as argument {2} (an empty string when there is no hint).
   Returns nullptr for TrendConclusion::NoConclusion.
   */
  const char *conclusion_txt_key( const TrendConclusion conclusion );

  /** Localization key for the confidence qualifier ("may", "probably", "strongly") used as
   argument {1} of the conclusion message.
   */
  const char *confidence_txt_key( const TrendConfidence confidence );

  /** The optional hint clause of the message: a localization key whose argument {1} is #value_txt
   (and {2} the localization key #zclass_key, when not null).  `key` is nullptr when there is no hint.
   */
  struct MessageHint
  {
    const char *key = nullptr;
    std::string value_txt;
    const char *zclass_key = nullptr;
  };

  MessageHint message_hint( const TrendResult &trend );

  /** Localization key of the sentence appended when TrendResult::possibleCascadeSumming. */
  const char *summing_caveat_txt_key();

  /** Plain (English) one-sentence conclusion text (plus the summing caveat), for non-GUI reports/logs
   (batch and inja templates) that are not localized.  Mirrors the GUI message.
   Returns an empty string for TrendConclusion::NoConclusion.
   */
  std::string conclusion_report_text( const TrendResult &trend );

  /** Areal density formatted for messages, e.g. "25 g/cm²" (two significant figures). */
  std::string areal_density_txt( const double areal_density );
}//namespace ShieldSourcePullTrend

#endif //ShieldSourcePullTrend_h
