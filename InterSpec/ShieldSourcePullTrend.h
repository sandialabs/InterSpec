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

#include <string>
#include <vector>
#include <utility>

// Forward declarations
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


/** Fits a low-order trend to the per-peak pulls ("chi" values) of an activity/shielding fit,
 and interprets the slope/curvature to diagnose likely model problems.

 Physics background: the pull for each peak is (observed-model)/uncertainty.  When the model
 shielding amount is wrong (and was not fit for), the attenuation-vs-energy mismatch produces
 a monotonic trend of the pulls with energy.  When the effective atomic number is wrong,
 photoelectric absorption (~Z^4.5/E^3) vs the nearly Z-independent Compton scattering means
 the *shape* of attenuation vs energy is wrong even after the fit compensates the overall
 amount - producing curvature (a low-energy "hook") in the pulls.

 Under a correct model each pull is ~N(0,1), so an unweighted polynomial least-squares fit
 gives coefficient t-statistics with an (approximately) known null distribution; thresholds
 on those t-statistics decide when a conclusion is warranted.

 For multiple nuclides a joint fit is performed: a separate intercept per nuclide (each
 nuclide's activity was independently fit, so its pulls center near zero), with slope and
 curvature shared across nuclides (all peaks see the same shielding-induced shape).
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
    /** The shielding's effective atomic number appears lower than reality.  The calibration
     study (scratch/chi_trend_study) confirmed the pull curvature (in log-energy) is concave-
     down when the model Z is too low - once poorly-converged fits are down-weighted via the
     residual-variance scaling.
     */
    AtomicNumberTooLow,
    /** The shielding's effective atomic number appears higher than reality (pull curvature
     concave-up).
     */
    AtomicNumberTooHigh,
    /** Both atomic number and areal density were fit, yet a trend remains - something
     else (distance, geometry, source age, detector response) may be off.
     */
    OtherModelInconsistency
  };//enum class TrendConclusion


  /** How strongly the test statistic supports the conclusion; controls message wording
   ("may be" / "probably" / "strongly suggests").  Tier boundaries come from the
   calibration study.
   */
  enum class TrendConfidence : int
  {
    /** A soft advisory ("possibly"): the observed trend is weak (t ~ 1.5-2.5) but the data
     demonstrably HAS the power to resolve the effect, so it is worth mentioning.  Only reachable
     for atomic-number calls, which are gated on the discrimination-power check.
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
  };//struct FitConfigSummary


  /** Builds a #FitConfigSummary from the shielding and source definitions of a fit. */
  FitConfigSummary summarize_fit_config(
                      const std::vector<ShieldingSourceFitCalc::ShieldingInfo> &shieldings,
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


  /** Result of the joint pull-trend fit. */
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

    std::vector<NuclideTrend> nuclideTrends;

    /** Atomic-number DISCRIMINATION POWER: the curvature t-statistic that a reference effective-Z
     error (a factor #sm_an_reference_z_factor change) WOULD produce for this exact peak set
     (energies + significances) and model shielding, if it existed.  It is a design property of
     the data, not of the observed pulls.

     Physics: above ~250 keV attenuation is Compton-dominated and nearly Z-independent in shape,
     so only low-energy, high-significance peaks carry an atomic-number signal.  A large value
     means the data could resolve a Z error; a small value (roughly < the power gate, calibrated
     to ~2.5) means it could not - so a null result is "inconclusive", not "atomic number
     consistent", and a curvature that DOES appear cannot be a Z signature (it is suppressed).

     Computed unit-robustly from the model per-peak optical depth (-ln of the shield attenuation
     fraction) and the mu(Z_ref)/mu(Z_model) ratio, so it does not depend on the fitted areal
     density.  0 when not computed (no shielding in the model, or fewer than a few points).
     */
    double anDiscriminationT = 0.0;

    /** True when #anDiscriminationT reaches the calibrated power gate, i.e. the data could actually
     resolve a reference-sized atomic-number error (calibration: detection of a 1.5x Z error
     becomes reliable at power ~2.5).  Atomic-number conclusions are only asserted when this is true.
     */
    bool anDiscriminable = false;

    TrendConclusion conclusion = TrendConclusion::NoConclusion;
    TrendConfidence confidence = TrendConfidence::Weak;
  };//struct TrendResult


  /** Fits the joint pull-vs-energy trend: per-nuclide intercepts plus shared slope
   (and shared quadratic term, when statistically determined), in centered/scaled
   `basis` coordinates.

   \param peak_comparisons Per-peak pull info, as from ShieldingSourceChi2Fcn::energy_chi_contributions
   \param peak_details Parallel per-peak details (for nuclide assignment and source fractions);
          must correspond index-for-index with `peak_comparisons` (extra entries in either are ignored)
   \param shieldings The model shieldings; only used when `basis` is TrendBasis::MuModel (may be empty otherwise)
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
   takes the qualifier text (see #confidence_txt_key) as argument {1}.
   Returns nullptr for TrendConclusion::NoConclusion.
   */
  const char *conclusion_txt_key( const TrendConclusion conclusion );

  /** Localization key for the confidence qualifier ("may", "probably", "strongly") used as
   argument {1} of the conclusion message.
   */
  const char *confidence_txt_key( const TrendConfidence confidence );

  /** Plain (English) one-sentence conclusion text, for non-GUI reports/logs (batch and inja
   templates) that are not localized.  Mirrors the GUI message (conclusion + confidence wording).
   Returns an empty string for TrendConclusion::NoConclusion.
   */
  std::string conclusion_report_text( const TrendResult &trend );
}//namespace ShieldSourcePullTrend

#endif //ShieldSourcePullTrend_h
