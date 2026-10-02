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

#ifndef FitPeaksCorpusScore_h
#define FitPeaksCorpusScore_h

/** Scoring library for the fit-peaks-for-nuclides corpus evaluation (fit_peaks_corpus_eval).

 Compares a set of fitted peaks against a manually fit reference set on the same spectrum, scoring
 peak presence, ROI grouping (share vs separate), continuum family, ROI extent, and fitted
 parameters with one canonical detection statistic.  Pure computation: no file output, so the
 terms can be unit tested on synthetic peak sets (test_fitPeaksCorpusScore.cpp).
 */

#include <map>
#include <string>
#include <vector>
#include <memory>

#include "InterSpec/PeakDef.h"
#include "InterSpec/RelActCalcAuto.h"

namespace SpecUtils
{
  class Measurement;
}

namespace FitPeaksCorpus
{

/** Continuum families the score compares.  CDF and non-CDF step variants map to the same family
 (the reference fits were meant to use the CDF forms throughout). */
enum class ContinuumFamily : int
{
  Linear = 0,     // NoOffset, Constant, Linear
  Poly2Plus,      // Quadratic, Cubic
  FlatStep,       // FlatStep, FlatStepCDF
  LinearStep,     // LinearStep, LinearStepCDF
  BiLinearStep,   // BiLinearStep, BiLinearStepCDF
  Other           // External
};

ContinuumFamily continuum_family( const PeakContinuum::OffsetType type );
const char *family_name( const ContinuumFamily family );
bool is_step_family( const ContinuumFamily family );

/** Fraction of a Gaussian's area within +/- `num_fwhm` FWHM of its mean; e.g. 0.9815 for 1 FWHM. */
double gaussian_fraction_within_num_fwhm( const double num_fwhm );

/** The canonical detection statistic z = S / sqrt(S + B) with S the Gaussian counts within
 +/-1 FWHM of the mean and B the peak's own continuum integrated over the same window (`roi_peaks`
 are needed for CDF step continua; `data` for data-defined steps).  Falls back to the gross data
 counts minus S when the peak has no continuum.  Returns 0 for non-Gaussian peaks. */
double detection_z( const PeakDef &peak,
                    const std::shared_ptr<const SpecUtils::Measurement> &data,
                    const std::vector<std::shared_ptr<const PeakDef>> &roi_peaks,
                    double *continuum_counts = nullptr );


/** One peak of a truth or fitted set, with the derived quantities the score uses. */
struct ScoredPeak
{
  std::shared_ptr<const PeakDef> peak;
  double energy = 0.0;
  double sigma = 0.0;
  double fwhm = 0.0;
  double amplitude = 0.0;
  double amplitude_uncert = 0.0;
  double continuum_counts = 0.0;   // B over +/-1 FWHM
  double z_det = 0.0;
  PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::NoOffset;
  ContinuumFamily family = ContinuumFamily::Other;
  size_t roi_index = 0;
  std::string source_name;         // nuclide/element symbol or reaction name; empty if none
  PeakDef::SourceGammaType gamma_type = PeakDef::NormalGamma;
  bool dont_care = false;          // truth peaks below the don't-care threshold
  int match = -1;                  // index into the other set, or -1
  int truth_cluster = -1;          // index into InjectTruth::merged_signal this peak was pooled into, or -1
  std::string verdict;             // see score_problem()
};//struct ScoredPeak


/** One ROI (shared PeakContinuum) of a truth or fitted set. */
struct ScoredRoi
{
  double lower = 0.0;
  double upper = 0.0;
  size_t first_channel = 0;
  size_t last_channel = 0;
  PeakContinuum::OffsetType continuum_type = PeakContinuum::OffsetType::NoOffset;
  ContinuumFamily family = ContinuumFamily::Other;
  std::vector<size_t> peaks;   // indices into PeakSet::peaks, sorted by energy
  size_t dominant = 0;         // index (into PeakSet::peaks) of the largest-amplitude peak
  double width_fwhm = 0.0;     // (upper - lower) / FWHM of the dominant peak
};//struct ScoredRoi


struct PeakSet
{
  std::vector<ScoredPeak> peaks;   // sorted by energy
  std::vector<ScoredRoi> rois;     // sorted by lower energy
};//struct PeakSet


/** Builds a scored peak set.  `truth_min_z` >= 0 marks peaks with z_det below it as don't-care
 (truth sets); pass a negative value for fitted sets.  Non-Gaussian peaks are skipped. */
PeakSet make_peak_set( const std::vector<std::shared_ptr<const PeakDef>> &peaks,
                       const std::shared_ptr<const SpecUtils::Measurement> &data,
                       const double truth_min_z );

/** Keeps only the peaks whose source name is in `sources` (used for single-NORM-nuclide scoring). */
std::vector<std::shared_ptr<const PeakDef>> filter_peaks_by_source(
  const std::vector<std::shared_ptr<const PeakDef>> &peaks,
  const std::vector<std::string> &sources );


/** One expected photopeak of the GADRAS-inject truth (a row of an "Expected ... Photopeaks"
 section of `<Source>_truth.csv`): the expectation value of the peak the detector would show for
 the gamma lines GADRAS merged into it. */
struct TruthPhotopeak
{
  double energy = 0.0;          // EffectiveEnergy (keV)
  double energy_lo = 0.0;       // energy range of the merged rows (equal to `energy` for one row)
  double energy_hi = 0.0;
  double area = 0.0;            // PeakArea: expected net counts
  double continuum_area = 0.0;  // ContinuumArea over [roi_lower, roi_upper]
  double fwhm = 0.0;            // EffectiveFWHM (keV)
  double roi_lower = 0.0, roi_upper = 0.0;
  double nsigma = 0.0;          // NSigmaOverBkg as written (PeakArea / sqrt(ContinuumArea))
  size_t num_lines = 0;         // NumGammaLines
  size_t num_rows = 1;          // rows merged into this entry (see merge_unresolved_photopeaks)
  bool signal = true;           // "Expected Signal Photopeaks" (else "Expected Background Photopeaks")

  /** The continuum expected within +/-1 FWHM of the energy (ContinuumArea scaled to that window). */
  double continuum_within_fwhm() const;
  /** The canonical detection statistic S/sqrt(S+B) with S = 0.9815 area, B = continuum_within_fwhm(). */
  double z_det() const;
};//struct TruthPhotopeak


/** The truth of one GADRAS-inject problem: the source's own expected photopeaks, the
 background's, and the resolution-merged signal list the scores use. */
struct InjectTruth
{
  std::string path;
  std::string source_name;                  // the "# Source:" header
  std::vector<TruthPhotopeak> signal;       // as listed, energy sorted
  std::vector<TruthPhotopeak> background;
  std::vector<TruthPhotopeak> merged_signal;      // rows closer than a mean FWHM merged
  std::vector<TruthPhotopeak> merged_background;
};//struct InjectTruth

/** Parses a `<Source>_truth.csv` (both photopeak sections).  Throws on unreadable files or a
 missing "Expected Signal Photopeaks" section; an empty section is not an error. */
InjectTruth parse_inject_truth_csv( const std::string &path );

/** Merges adjacent rows closer than the mean of their FWHMs (area-summed, area-weighted energy,
 widest ROI), since two such lines are one peak in the data and in any fit. */
std::vector<TruthPhotopeak> merge_unresolved_photopeaks( const std::vector<TruthPhotopeak> &rows );

/** The inject problem name for a hand-fit corpus id: "q115-Lu177m-Unsh" -> "Lu177m_Unsh",
 "q098-K40-Sh-Point" -> "K40_Sh-Point" (the first separator only; the inject names keep the rest). */
std::string inject_name_for_problem_id( const std::string &problem_id );


/** Objective weights; every value is printed into run outputs so results are self-describing. */
struct ScoreWeights
{
  double min_scored_energy = 0.0;    // peaks below this energy are ignored entirely (both sets)
  double truth_min_z = 1.5;          // truth peaks below this are don't-care
  double match_num_fwhm = 1.0;       // match window in truth FWHM ...
  double match_min_kev = 0.5;        // ... but at least this many keV

  double weak_z = 3.0;               // truth classes: weak [truth_min_z, weak_z), moderate [weak_z, strong_z), strong
  double strong_z = 8.0;
  double missed_weak = 0.5;
  double missed_moderate = 1.5;
  double missed_strong = 4.0;

  double ghost_z = 1.0;              // fitted extras: ghost [0, ghost_z), weak [ghost_z, weak_z), significant
  double extra_ghost = 0.1;
  double extra_weak = 0.5;
  double extra_significant = 1.5;

  double share_max_sep_fwhm = 8.0;   // adjacent matched truth peaks closer than this form a pair
  double share_disagree = 1.5;

  double family_step_disagree = 1.0; // step vs non-step, or different step kinds
  double family_poly_disagree = 0.5; // Linear vs Poly2Plus

  double extent_per_fwhm = 0.5;      // per ROI side, beyond the dead band, capped
  double extent_deadband_fwhm = 0.5;
  double extent_cap_fwhm = 3.0;

  double area_pull_weight = 0.1;
  double area_rel_floor = 0.05;
  double area_pull_cap = 5.0;

  double mean_offset_weight = 0.5;
  double mean_offset_deadband_fwhm = 0.1;
  double mean_offset_cap_fwhm = 1.0;

  double failure = 20.0;
  double nondeterminism = 5.0;

  // Inject-truth terms (only when an InjectTruth is supplied to score_problem):
  double legit_truth_min_z = 1.0;    // an unmatched fitted peak on a signal photopeak at least this significant is "legit", not an extra
  /** A truth photopeak the FOREGROUND itself does not show is not a miss.  The inject truth is
   computed from the source term, not from the recorded spectrum, so it lists peaks a given detector
   never sees - the Fulcrum 40h has no useful response below ~70 keV, yet its truth files claim
   tens of thousands of counts at 20-40 keV.  Charging the fitter for not finding a peak that is
   absent from the data measures the truth file, not the fitter.  <= 0 disables the test. */
  double miss_requires_data_z = 3.0;
  /** Grade a missed truth peak by what the SPECTRUM shows, not only by the area the truth file
   claims, taking the smaller of the two.  See `miss_requires_data_z` for why the claims cannot be
   trusted on their own. */
  bool class_miss_by_data_z = true;
  /** Shielding fluoresces, and the inject truth lists those x-rays among the source's own
   photopeaks: four of Xe133_Sh's strong "missed" peaks are the complete Pb K series (72.80, 74.97,
   84.94, 87.36 keV, matched to better than 0.15 keV).  No requested nuclide can produce them - the
   emitter is the shield - so charging the fitter for missing them measures the geometry, not the
   fit.  Excused only when no requested source has a line there.  <= 0 disables. */
  double shield_fluorescence_num_fwhm = 1.0;
  /** Annihilation is real and detector-wide, but the GADRAS photopeak lists never enumerate it, so a
   fitted 511 keV peak would always score as an extra.  Free, like "xray". */
  bool free_annihilation_peak = true;
  double truth_area_weight = 0.0;    // weight of the fitted-vs-truth area pull per truth cluster (0 = report only)
  double truth_area_rel_floor = 0.03;// relative systematic floor on the truth area
  double truth_area_cap = 5.0;       // |pull| cap
  double truth_area_min_z = 3.0;     // clusters below this z are graded only when a peak was fit on them
  double legit_max_pull = 4.0;       // a "legit" extra may not exceed the truth area by more than this many sigma
  double merged_max_sep_fwhm = 1.0;  // a reference peak is only excused as "merged" when a MATCHED reference peak sits this close

  std::string to_string( const std::string &separator = " " ) const;
};//struct ScoreWeights


/** An adjacent matched truth-peak pair and the share/separate verdicts on each side. */
struct PairRecord
{
  size_t truth_a = 0, truth_b = 0;
  double separation_fwhm = 0.0;
  bool truth_share = false;
  bool fit_share = false;
};

/** A truth ROI compared to the fitted ROI containing the match of its dominant peak. */
struct RoiRecord
{
  size_t truth_roi = 0;
  int fit_roi = -1;
  ContinuumFamily truth_family = ContinuumFamily::Other;
  ContinuumFamily fit_family = ContinuumFamily::Other;
  double d_lower_fwhm = 0.0;   // (fit.lower - truth.lower) / FWHM
  double d_upper_fwhm = 0.0;
  double cost_family = 0.0;
  double cost_extent = 0.0;
};


/** One resolution-merged truth photopeak and the areas the fitted and the reference sets put on it. */
struct TruthAreaRecord
{
  size_t cluster = 0;           // index into InjectTruth::merged_signal
  double energy = 0.0;
  double area = 0.0;            // expected net counts
  double z = 0.0;               // TruthPhotopeak::z_det()
  double sigma = 0.0;           // sqrt( area + continuum within +/-1 FWHM + (rel_floor*area)^2 )
  size_t fit_n = 0, ref_n = 0;  // peaks pooled into the cluster
  double fit_area = 0.0, ref_area = 0.0;
  double fit_pull = 0.0, ref_pull = 0.0;   // (pooled area - truth area) / sigma; NaN when not graded
  bool fit_graded = false, ref_graded = false;
  bool xray_affected = false;   // x-rays contribute here: not graded on area (yields are unreliable)
};


struct ProblemScore
{
  size_t n_truth_scored = 0;
  size_t n_truth_dontcare = 0;
  size_t n_fitted = 0;
  size_t n_matched = 0;          // scored truth peaks with a match
  size_t missed_weak = 0, missed_moderate = 0, missed_strong = 0;
  size_t extra_ghost = 0, extra_weak = 0, extra_significant = 0;
  size_t neutral = 0;            // fitted peaks matching a don't-care truth peak
  size_t pairs = 0, share_disagree = 0;
  size_t rois_compared = 0, family_disagree = 0;
  size_t extent_sides = 0, extent_gt1fwhm_sides = 0;

  // Inject-truth terms (zero unless score_problem was given an InjectTruth)
  size_t truth_on_bkg_line = 0;  // reference peaks on a background photopeak with no signal photopeak (don't-care)
  size_t truth_shield_xray = 0;       // truth photopeaks that are shielding fluorescence, not source lines
  size_t truth_absent_from_data = 0;  // truth photopeaks the foreground itself does not show (don't-care)
  size_t fitted_annihilation = 0;     // fitted 511 keV peaks the photopeak lists never enumerate (free)
  size_t truth_merged_away = 0;  // reference peaks the truth merges into a photopeak the fit did report
  size_t extra_legit = 0;        // unmatched fitted peaks that sit on a signal photopeak (not penalised)
  size_t extra_xray = 0;         // ... or on a characteristic x-ray of a requested source (not penalised)
  size_t extra_bkg_line = 0;     // unmatched fitted peaks that sit on a background photopeak (penalised as extras)
  size_t truth_clusters = 0;     // merged signal photopeaks
  size_t truth_xray_clusters = 0;// ... of which x-rays contribute to, so not graded on area
  size_t truth_fit_graded = 0, truth_ref_graded = 0;
  size_t truth_fit_bad = 0, truth_ref_bad = 0;          // |pull| > 3
  double truth_fit_abs_pull = 0.0, truth_ref_abs_pull = 0.0;   // sums of capped |pull|

  double cost_missed = 0.0, cost_extra = 0.0, cost_share = 0.0, cost_family = 0.0;
  double cost_extent = 0.0, cost_area = 0.0, cost_mean = 0.0, cost_failure = 0.0;
  double cost_truth_area = 0.0;

  std::vector<PairRecord> pair_records;
  std::vector<RoiRecord> roi_records;
  std::vector<TruthAreaRecord> truth_area_records;

  double raw_cost() const;
  double norm_cost() const;   // raw_cost / max(1, n_truth_scored)
};//struct ProblemScore


/** Scores `fitted` against `truth`, filling the match/verdict fields of both sets.
 Truth verdicts: "matched", "missed", "dontcare" (unmatched don't-care), "matched_dontcare", and
 with an `inject` truth "dontcare_bkg" (a reference peak on a background photopeak where the source
 has no photopeak: the reference label is wrong there, so the peak is neither missed nor matched) and
 "dontcare_merged" (an unmatched reference peak that is unresolvable - within `merged_max_sep_fwhm`
 - from a MATCHED reference peak on the same truth photopeak: the truth says there is one detectable
 peak there and the fit reported it, so the reference having split it is not a miss.  Sharing a truth
 cluster is not enough on its own: the truth's own merge is single-linkage and can chain across
 several FWHM, and a cluster that wide holds peaks a good fit should resolve).
 Fitted verdicts: "matched", "neutral", "ghost", "extra_weak", "extra", and with an `inject` truth
 also "legit" (an unmatched peak on a signal photopeak of at least `legit_truth_min_z`; costs
 nothing - the reference simply did not fit it) and "extra_bkg" (on a background photopeak; a real
 peak wrongly attributed to the source, costed like an extra).  Also "xray": an unmatched peak on a characteristic x-ray of a requested source is real, and
 charging for it would tune the fitter away from finding source x-rays.  The inject truth grades the
 pooled fitted and reference areas on each merged signal photopeak (`truth_area_records`). */
ProblemScore score_problem( PeakSet &truth, PeakSet &fitted, const ScoreWeights &weights,
                            const InjectTruth *inject = nullptr,
                            const std::vector<double> *source_xrays = nullptr,
                            const std::shared_ptr<const SpecUtils::Measurement> &foreground = nullptr,
                            const std::vector<double> *source_lines = nullptr );

/** Every gamma and x-ray energy the requested sources emit above `min_intensity` of their strongest,
 used to tell a peak the sources CAN explain from one only the surroundings can. */
std::vector<double> source_line_energies( const std::vector<RelActCalcAuto::SrcVariant> &sources,
                                          const double min_intensity = 1.0e-4 );


/** One reference problem: the spectra, the manually fit reference peaks, and the sources to fit. */
struct CorpusProblem
{
  std::string id;
  std::string path;
  std::shared_ptr<const SpecUtils::Measurement> foreground;
  std::shared_ptr<const SpecUtils::Measurement> background;   // may be null
  int foreground_sample = -1;
  std::vector<std::shared_ptr<const PeakDef>> truth_peaks;    // deep copies with private continua
  std::vector<std::string> requested_source_names;
  std::vector<RelActCalcAuto::SrcVariant> sources;
  std::shared_ptr<const InjectTruth> inject_truth;   // GADRAS truth for the same spectrum, if attached
  std::vector<std::string> notes;   // parse warnings, special-case handling
};//struct CorpusProblem

/** Attaches the GADRAS-inject truth for a hand-fit corpus problem from `inject_dir` (a
 `<Detector>/<Location>/<time>_seconds` directory): parses `<name>_truth.csv` and checks that the
 problem's foreground is channel-for-channel the `<name>.pcf` record 0 (a note records a mismatch,
 and the truth is still attached).  Returns false when no truth file exists for the id. */
bool attach_inject_truth( CorpusProblem &problem, const std::string &inject_dir );

/** Problem id from a file path: the file name without extension ("q054-Cs137-Unsh"). */
std::string problem_id_from_path( const std::string &path );

/** Splits "Pu, Pu238,Pu239" into trimmed, non-empty names. */
std::vector<std::string> split_source_list( const std::string &text );

/** Loads an InterSpec-written N42 (foreground + optional background, embedded user peaks).
 Requested sources come from a foreground "Requested sources:" remark, else `manifest_sources`
 (id -> comma list), else the foreground title.  `InterSpec::setStaticDataDirectory` must have
 been called first so nuclide assignments decode.  Throws on unreadable files. */
CorpusProblem load_problem( const std::string &path,
                            const std::map<std::string,std::string> &manifest_sources );

/** Characteristic x-ray energies (keV) of the requested sources: the x-rays emitted by a nuclide's
 own decay chain (SandiaDecay), and, for an element source, that element's fluorescence x-rays.
 A fitted peak on one of these is a real feature of the source even though neither the hand fits nor
 the GADRAS "expected signal photopeaks" list x-rays - the Cd109 truth file, for instance, contains
 only the 88 keV gamma while the spectrum plainly shows the Ag K x-rays at 22 and 25 keV.  Only
 entries carrying at least `min_intensity` of the strongest x-ray are returned. */
std::vector<double> source_xray_energies( const std::vector<RelActCalcAuto::SrcVariant> &sources,
                                          const double min_intensity = 0.01 );

/** Loads one problem of the GADRAS-inject dataset (peak_fit_accuracy_inject_compact).
 `truth_csv_path` is a `<Source>_truth.csv`; the spectra come from the sibling `<Source>.pcf`
 (record 0 = Poisson-varied source + background, record 1 = a background of the same duration).
 Truth peaks are the "Expected Signal Photopeaks" rows: one Gaussian per effective photopeak
 (EffectiveEnergy, PeakArea, EffectiveFWHM) on a flat continuum carrying ContinuumArea over
 [ROILower, ROIUpper].  The ROI structure is synthetic, so score these with the structure weights
 zeroed (`--no-structure`).  The requested source is the file name's first token (Am241Li -> Am241,
 Tl201wTl202 -> Tl201,Tl202, Uore -> U238,U235).  Throws on unreadable files or unknown sources. */
CorpusProblem load_inject_problem( const std::string &truth_csv_path );


/** The result of fitting one problem under one variant, everything the reports need. */
struct ProblemResult
{
  std::string id;          // problem id, plus a variant suffix when applicable
  /** Full path of the spectrum this problem was loaded from, so the gallery can name it and the
   reader can open the file without hunting through the corpus directory. */
  std::string spectrum_path;
  /** ROIs the planner produced, before the solve; see PeakFitResult::planned_rois. */
  std::vector<std::pair<double,double>> planned_rois;
  std::string mode;
  std::string sources;     // comma-separated requested sources
  std::string status;
  std::string error;
  std::vector<std::string> warnings;
  bool mechanical_failure = false;
  bool nondeterministic = false;
  size_t dev_check_failures = 0;
  double wall_seconds = 0.0;
  PeakSet truth;
  PeakSet fitted;
  ProblemScore score;
  std::shared_ptr<const SpecUtils::Measurement> foreground;
  std::shared_ptr<const SpecUtils::Measurement> background;
  std::vector<std::shared_ptr<const PeakDef>> fitted_peaks;   // the set that was scored
  std::vector<std::shared_ptr<const PeakDef>> truth_peaks;
  PeakSet autosearch;       // the automated-search peaks the fit was given (diagnostics only)
  std::shared_ptr<const InjectTruth> inject_truth;   // as scored against (may be null)
  std::vector<std::string> roi_plan_trace;   // planner decisions (empty unless use_roi_plan)
  /** True when `truth` came from a GADRAS inject truth file rather than from hand-fit reference
   peaks.  The two are scored the same way but must not be *labelled* the same: calling a computed
   photopeak list "reference peaks" invites reading a truth-file artifact as the user's own fit. */
  bool truth_is_inject = false;
};//struct ProblemResult

/** A fingerprint of the fitted peaks (energies, areas, ROI bounds) for determinism checks. */
std::string fitted_fingerprint( const PeakSet &fitted );

}//namespace FitPeaksCorpus

#endif //FitPeaksCorpusScore_h
