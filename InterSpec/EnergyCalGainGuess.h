#ifndef EnergyCalGainGuess_h
#define EnergyCalGainGuess_h
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

#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include <cstdint>

#include <Wt/WDialog.h>
#include <Wt/WContainerWidget.h>

// Forward Declarations
class InterSpec;
class AuxWindow;
class EnergyCalTool;
class D3SpectrumDisplayDiv;
struct MeasToApplyCoefChangeTo;

namespace Wt
{
  class WText;
  class WTable;
  class WTableRow;
  class WLineEdit;
  class WPushButton;
  class WSuggestionPopup;
}

namespace SpecUtils
{
  class Measurement;
  struct EnergyCalibration;
}

namespace SandiaDecay
{
  struct Nuclide;
  struct Element;
}


/** GUI-free computation to guess the energy calibration (linear offset + gain) of a spectrum
 that has no usable calibration, so unit tests (target/testing/test_EnergyCalGainGuess.cpp) can
 exercise each stage without a WApplication.

 The overall approach automates what a spectroscopist does by hand with an uncalibrated
 spectrum:
   1) Find candidate peaks purely in channel space (an identity calibration is installed on a
      copy of the spectrum, so the second-derivative peak finder returns channel units).
   2) Generate many (offset=0, gain) hypotheses:
        - Pairs of strong candidate peaks whose channel *ratio* matches the energy ratio of two
          known gamma/x-ray lines (a gain-invariant fingerprint when offset ~ 0), which covers
          the familiar-peak-spacing, fluorescence Ka/Kb, and highest-peak heuristics at once.
        - Single anchors: highest-energy peak as 2614.5/1460.8 keV, biggest peak as
          K-40/Cs-137/etc, anomalously-wide peaks as 511 keV.
        - The spectrums right edge at common endpoint energies (400 keV ... 10 MeV).
        - Common gains (e.g. 1 keV/channel), and a coarse full-range gain scan as a backstop.
   3) Score every hypothesis with one common metric: how much of the significant peak
      significance is explained by known lines - discounted by the chance of accidental
      matches - plus nuclide self-consistency, FWHM-vs-energy plausibility (a 511 keV match is
      *expected* to be wide), escape/backscatter credit, and penalties for unexplained peaks.
   4) The best-scoring hypotheses are refined by a linear (offset,gain) least-squares fit to
      their peak-line assignments, re-scored, de-duplicated, and returned ranked.

 Backscatter *shape* template matching was considered but not implemented; backscatter peaks of
 matched lines are instead credited in the scoring, which captures most of the value.
 */
namespace GainGuessCalc
{
  /** A candidate photopeak found in channel space (no trustworthy energy calibration). */
  struct CandidatePeak
  {
    /** Peak mean, in (fractional) channels. */
    double channel = 0.0;

    /** Gaussian sigma, in channels. */
    double sigma_ch = 0.0;

    /** Peak area, in counts, and its (approximate, Poisson-ish) uncertainty. */
    double area = 0.0;
    double area_uncert = 0.0;

    /** True for the strongest peaks - the ones a hypothesis is expected to explain. */
    bool significant = false;

    /** True if this peak is notably wider than the width-vs-channel trend of the other peaks
     (only evaluated for high-resolution spectra) - a natural 511 keV / x-ray suspect.
     */
    bool width_outlier = false;
  };//struct CandidatePeak


  /** One known emission line a peak can be matched to. */
  struct LineSource
  {
    enum class Type : int
    {
      NuclideGamma,
      Annihilation,
      FluorXray,     //fluorescence x-ray of an element (shielding/self-attenuating source)
      SingleEscape,  //only produced in the scorer, from a matched primary
      DoubleEscape,  // "
      Backscatter    // "
    };//enum class Type

    Type type = Type::NuclideGamma;

    /** Line energy, keV. */
    double energy = 0.0;

    /** Emission intensity relative to the strongest line of the same parent, in (0,1]. */
    float rel_intensity = 0.0f;

    /** Parent nuclide (nullptr for Annihilation/FluorXray). */
    const SandiaDecay::Nuclide *nuc = nullptr;

    /** Fluorescing element (FluorXray only). */
    const SandiaDecay::Element *el = nullptr;

    /** True if this line came from a user-supplied candidate nuclide (gets a score bonus). */
    bool user_supplied = false;
  };//struct LineSource


  /** The library of lines hypotheses are scored against, plus the precomputed pairwise
   energy-ratio table used by the pair-matching hypothesis generator.
   */
  struct LineLibrary
  {
    /** All lines, sorted by energy. */
    std::vector<LineSource> lines;

    /** (E_low/E_high, index of low line, index of high line), sorted by ratio; only pairs of
     reasonably strong lines are included.
     */
    struct RatioEntry
    {
      float ratio;
      uint32_t low_index;
      uint32_t high_index;
    };
    std::vector<RatioEntry> ratios;
  };//struct LineLibrary


  /** A candidate-peak <-> known-line assignment of a scored hypothesis. */
  struct PeakAssignment
  {
    /** Index into the CandidatePeak vector. */
    size_t peak_index = 0;

    LineSource line;

    /** (predicted energy of peak) - (line energy), keV. */
    double delta_kev = 0.0;
  };//struct PeakAssignment


  /** Which heuristic produced the seed a #GuessedCal grew from (for display). */
  enum class SeedStrategy : int
  {
    PeakPair,     //channel ratio of two peaks matched an energy ratio of two known lines
    SingleAnchor, //a single prominent peak anchored to a specific line
    Wide511,      //an anomalously-wide peak taken as the 511 keV annihilation line
    Endpoint,     //the right edge of the data placed at a common endpoint energy
    CommonGain,   //a common keV/channel value
    GainScan      //the coarse full-range gain scan backstop
  };//enum class SeedStrategy


  /** One ranked result: a plausible linear energy calibration. */
  struct GuessedCal
  {
    /** E = offset + gain*channel; i.e. the two polynomial coefficients. */
    double offset = 0.0;
    double gain = 0.0;

    /** The heuristic that seeded this hypothesis, and (strategy-dependent) the line/endpoint
     energies involved - e.g. the two line energies of a PeakPair match.
     */
    SeedStrategy strategy = SeedStrategy::GainScan;
    double strategy_energy1 = 0.0;
    double strategy_energy2 = 0.0;

    /** The score - higher is better; only meaningful relative to other hypotheses for the
     same spectrum.
     */
    double score = 0.0;

    /** Fraction (0..1) of the significant peaks weight explained by known lines - the
     user-facing confidence-ish number.
     */
    double explained_frac = 0.0;

    std::vector<PeakAssignment> assignments;

    /** Implied parent nuclides, ordered by their contribution to the score. */
    std::vector<const SandiaDecay::Nuclide *> nuclides;

    /** Human-readable seed strategy, e.g. "peak pair 1460.8/2614.5 keV" or "wide peak as 511".
     Not translated - for logging/debugging; the GUI maps seed strategies to tr() keys itself.
     */
    std::string provenance;
  };//struct GuessedCal


  struct GuessOptions
  {
    /** Nuclides the user thinks may be present; their lines are added to the library (if not
     already there) with a scoring bonus.
     */
    std::vector<const SandiaDecay::Nuclide *> user_nuclides;

    /** Maximum number of ranked results to return. */
    size_t max_results = 10;

    /** Optional cancellation token - when non-null and set, computation stops at the next
     cooperative check and an empty result is returned.
     */
    std::shared_ptr<const std::atomic<bool>> cancel;
  };//struct GuessOptions


  /** A #GuessedCal whose #GuessedCal::explained_frac is below this is presented as
   low-confidence ("no confident guess") by the GUI.
   */
  constexpr double sm_min_confident_explained_frac = 0.5;


  /** Finds candidate peaks purely in channel space: installs an identity (channel-unit)
   calibration on a copy of the spectrum and runs the second-derivative candidate search, so
   the results do not depend on whatever (possibly nonsense) calibration the file carried.

   Keeps the ~30 most significant candidates; the strongest ~15 get #CandidatePeak::significant.
   Also classifies the spectrum as high (HPGe) or low/medium resolution from the relative
   candidate widths + channel count, and (for high-res) flags width outliers.

   @param meas The spectrum; must have at least 64 channels.
   @param is_high_res [out] whether the spectrum looks like HPGe data.
   @returns candidates sorted by channel; empty if the spectrum is unusable/featureless.

   Throws std::exception if meas is null or has too few channels.
   */
  std::vector<CandidatePeak>
  find_candidate_peaks_channelspace( const std::shared_ptr<const SpecUtils::Measurement> &meas,
                                     bool &is_high_res );


  /** Builds the line library: a curated common-nuclide list (NORM, medical/industrial sources,
   and SNM: K-40, Th-232 chain, U-238/U-235 chains, Np-237, Pu-238/239/240/241, Cs-137, Co-60,
   etc - each decayed to a representative age via SandiaDecay), fluorescence x-rays of
   W/Pb/Bi/U/Th, the 511 keV annihilation line, and any user-supplied nuclides.

   Throws std::exception if the decay database is not initialized.
   */
  LineLibrary build_line_library( const std::vector<const SandiaDecay::Nuclide *> &user_nuclides );


  /** Scores the linear-calibration hypothesis E = offset + gain*channel against the candidate
   peaks; see the namespace comment for the ingredients.  Higher is better; a hypothesis that
   explains nothing scores <= 0.  Returns -999 for hypotheses failing sanity gates (implied
   energy span, absurd offset).

   @param assignments If non-null, receives the peak-line assignments made.
   @param explained_frac If non-null, receives the explained fraction of significant peak weight.
   */
  double score_hypothesis( const double offset, const double gain,
                           const std::vector<CandidatePeak> &peaks,
                           const size_t nchannels,
                           const LineLibrary &library,
                           const bool is_high_res,
                           std::vector<PeakAssignment> *assignments = nullptr,
                           double *explained_frac = nullptr );


  /** Top-level driver: candidate search, hypothesis generation, scoring, refinement, ranking.

   @returns Up to GuessOptions::max_results ranked calibrations - possibly none, if the spectrum
            is featureless or nothing matched.  Results are returned even when low-confidence
            (see #sm_min_confident_explained_frac); the caller decides how to present them.

   Safe to call on a worker thread (pure computation; only touches the passed-in data and the
   decay database singleton).  Never throws; on error (null/short spectrum, uninitialized decay
   database, cancellation) returns an empty vector.
   */
  std::vector<GuessedCal>
  guess_energy_cal( const std::shared_ptr<const SpecUtils::Measurement> &meas,
                    const GuessOptions &options );
}//namespace GainGuessCalc


/** Content widget of the "Guess Gain..." dialog opened from the Energy Calibration tabs
 "More Actions" column; created by EnergyCalAddActionsWindow, with buttons in the AuxWindow
 footer (analogous to EnergyCalGainMatch).

 Runs GainGuessCalc::guess_energy_cal on the displayed foreground (on a worker thread), lists
 the ranked candidate calibrations, previews the selected one on an embedded spectrum chart
 (with the implied nuclides in the row), and lets the user add candidate-nuclide hints and
 re-run.  On "Use", the selected calibration is applied through
 EnergyCalTool::applyCalChange() (undoable), and the Reference Photopeaks display is set to the
 matched nuclide(s).
 */
class EnergyCalGainGuess : public Wt::WContainerWidget
{
public:
  EnergyCalGainGuess( std::shared_ptr<std::vector<MeasToApplyCoefChangeTo>> measToChange,
                      EnergyCalTool *cal, AuxWindow *parent );
  virtual ~EnergyCalGainGuess();

  void handleFinish( Wt::DialogCode result );

  /** Session-thread continuation of #startCalc; public so the worker post-back (which
   re-looks-up this widget by DOM id) can call it.  `generation` guards against stale results.
   */
  void onCalcDone( int generation, const std::vector<GainGuessCalc::GuessedCal> &results );

protected:
  void startCalc();                 //snapshots inputs and runs guess_energy_cal on a worker thread
  void handleAddNuclide();          //adds the typed nuclide as a chip + re-runs
  void removeNuclide( const SandiaDecay::Nuclide *nuc );
  void rebuildResultTable();
  void handleRowSelected( const int row_index );
  void updatePreview();
  void showPreviewRefLines( const GainGuessCalc::GuessedCal &guess );  //ref lines on the preview chart
  void applyChanges();

  InterSpec *m_interspec;
  EnergyCalTool *m_calibrator;
  AuxWindow *m_parent;
  std::shared_ptr<std::vector<MeasToApplyCoefChangeTo>> m_measToChange;

  /** The foreground spectrum snapshot the current results were computed for. */
  std::shared_ptr<const SpecUtils::Measurement> m_foreground;

  std::vector<const SandiaDecay::Nuclide *> m_userNuclides;
  std::vector<GainGuessCalc::GuessedCal> m_results;
  int m_selectedRow;                // -1 == none
  int m_generation;                 //bumped per startCalc; the post-back discards stale results
  bool m_calcRunning;

  /** Cancellation token for the in-flight worker computation (a fresh one per startCalc). */
  std::shared_ptr<std::atomic<bool>> m_cancelCalc;

  D3SpectrumDisplayDiv *m_chart;
  Wt::WLineEdit *m_nuclideEdit;
  Wt::WSuggestionPopup *m_nuclideSuggest;
  Wt::WContainerWidget *m_nuclideChips;
  Wt::WTable *m_resultTable;
  std::vector<Wt::WTableRow *> m_resultRows;  //data rows, parallel to m_results
  Wt::WText *m_status;
  Wt::WPushButton *m_cancel;
  Wt::WPushButton *m_use;
};//class EnergyCalGainGuess

#endif //EnergyCalGainGuess_h
