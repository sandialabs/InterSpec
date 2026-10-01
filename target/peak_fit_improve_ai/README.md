# peak_fit_improve_ai: tools for improving `fit_peaks_for_nuclides`

The harness, scripts and review material used, in AI-assisted sessions, to measure, look at and improve
`FitPeaksForNuclides::fit_peaks_for_nuclides` against peak-fit corpora kept outside the repository.

| Path | What it is |
|---|---|
| `harness/` | `fit_peaks_corpus_eval`, which fits a corpus and scores it against hand fits and/or GADRAS-inject truth, plus its scorer (`FitPeaksCorpusScore.*`), outputs and HTML gallery (`FitPeaksCorpusReport.*`), and the scorer's unit test (ctest `TFitPeaksCorpusScore`). |
| `tools/` | Shell scripts to snapshot binaries and run corpora, and Python scripts to compare runs, diagnose and render fits.  Tables below. |
| `review/` | Instructions for visual review by agents: `RUBRIC.md` (one fit; written for NaI), `RUBRIC_BOTH.md` (new vs old) and `calibration/` example images. |
| `harness/peak_fit_objective_eval.cpp` | Compares PeakFitLM's fit objectives (chi2, the default sparse-data likelihood refit, the likelihood forced on every ROI) on Poisson replicas with known truth - see the section below. |
| `CMakeLists.txt` | The targets; added from `target/testing/CMakeLists.txt`. |
| `*PROMPT*.md` | Session prompts (tracked, although the repository's `*prompt*.md` ignore rule matches them, so a new one needs `git add -f`): `HPGE_PROMPT.md` (the same effort for HPGe), `LOWRES_FOLLOWON_PROMPT.md` (continue the scintillator work), `INVESTIGATION_PROMPT.md` (the generic template for any other detector class or dataset). |

## Build and setup

- **Build** as part of the test tree (Release with `PERFORM_DEVELOPER_CHECKS` on):
  `cd target/testing/build_ninja && ninja fit_peaks_corpus_eval test_fitPeaksCorpusScore`.
  The executables land at the top of that build tree.  See `fit_peaks_corpus_eval --help`.
- **Python** for the image scripts needs numpy and matplotlib, which the system `python3` lacks:
  `python3 -m venv ~/fit_peaks_work/venv && ~/fit_peaks_work/venv/bin/pip install numpy matplotlib`,
  then `export FPR_PY=~/fit_peaks_work/venv/bin/python`.
- **Settings.** All shell scripts source `tools/env.sh`; override its variables from the environment:
  `FPR_WORK` (work area, default `~/fit_peaks_work`: `bins/`, `runs/`, `ov/`, `triage/`, `galleries/`),
  `FPR_INJ` (GADRAS-inject corpus root), `FPR_SITE` / `FPR_DWELL` (inject background site and dwell,
  default Livermore / 300), `FPR_MANUAL` (R500 NaI hand fits), `FPR_DATA`, `FPR_PY`, `FPR_THREADS`.

## Corpora

- **GADRAS-inject** (`$FPR_INJ/<detector>/<site>/<dwell>_seconds/<Source>.pcf` + `<Source>_truth.csv`):
  223 source spectra per detector for 35 detectors, Livermore/Baltimore/Denver backgrounds, 30/300/1800 s.
  Scored for presence and area against the truth (`--corpus-format=inject --no-structure`).  Used so far:
  IdentiFINDER-R500-NaI, IdentiFINDER-NGH, SAM-Eagle-NaI-3x3, Radseeker-LaBr3, Kromek-GR1-CZT (det-types
  `Low`/`LaBr`/`CZT`) and Detective-X, Detective-EX, HPGe_Planar_50%, Fulcrum40h (`High`).
- **Hand fits**: the 168 Detective-X 300 s HPGe spectra are the harness's default corpus (reference ROIs,
  sharing and continua scored too, with the inject truth attached); `$FPR_MANUAL/IdentiFINDER-R500-NaI_300_seconds`
  holds 17 R500 NaI hand fits (`run_hand.sh`).
- **Always** pass `--det-type` for scintillators (the harness defaults to `High`), `--background=file`,
  and `--weight=min_scored_energy=20`; the scripts do.

## Running the fitter
| Script | What it does |
|---|---|
| `snapshot.sh TAG` | Freezes the current build as `bins/eval_TAG`, with its sources and diff in `bins/src_TAG/`. The build tree is shared, so measure only snapshots. |
| `run_corpus.sh BIN OUT DETECTOR DET_TYPE [TIMEOUT] [args]` | One full inject-corpus run. `DET_TYPE` is `Low`, `LaBr`, `CZT` or `High`. |
| `quick.sh BIN OUT DETECTOR DET_TYPE a,b,c [args]` | Runs a subset of problems and renders images. Set `FPR_CMP=<run dir>` to draw the old fit in blue. |
| `full.sh BIN TAG [args]` | All five scintillator sets, Detective-X and the HPGe hand fits, in two parallel streams (1.5-2 h). Extra args (e.g. `--set field=value`) go to every run. With `FPR_GUARD_REF=<tag>` it checks bit-identity of the runs in `FPR_GUARD_SETS` (default `detx manual`). |
| `full_hpge.sh BIN TAG [args]` | The HPGe sets: Detective-X at 30/300/1800 s, Detective-EX, HPGe_Planar_50%, Falcon 5000, Fulcrum40h and the hand fits. |
| `guard.sh BIN TAG REF [args]` | The HPGe guard alone (Detective-X + hand fits) with per-peak bit-identity against `REF`. |
| `guard_lowres.sh BIN TAG REF [args]` | The scintillator guard (the five sets) with per-peak bit-identity against `REF`: run it after HPGe changes. |
| `gallery.sh BIN NAME [DET:TYPE ... \| manual]` | The harness's HTML galleries into `galleries/NAME/<detector>/gallery.html`; `manual` = the HPGe hand fits. |
| `run_hand.sh BIN OUT` | The 17 hand-fit R500 NaI spectra, scored against the hand fits. |

## Comparing runs
| Tool | What it does |
|---|---|
| `cmp_runs.py RUN_A RUN_B` | Truth lines found, split into strong (z>=8), moderate and weak, plus every line gained or lost. The main A/B tool. Compare these counts, not raw costs, across scorer changes. |
| `run_table.py PREFIX...` | Strong/moderate/weak found and significant extra peaks per detector of `full.sh` runs (`PREFIX_r500` ... `PREFIX_czt`), side by side. |
| `extras_cmp.py RUN_A RUN_B` | Significant extra peaks in each run, and the ones B adds, strongest first. Read it whenever a change finds more lines. |
| `miss_why.py RUN [minz] [maxz]` | Why each missed truth line was missed: covered by the last plan, or rejected (with the planner's reason). |
| `invariants.py RUN...` | Things that must never happen: fitted peaks outside their own ROI, overlapping ROIs, ghost peaks (z<1), timeouts. |
| `bitident.py RUN_A RUN_B` | Checks that every fitted peak is identical. |
| `compare_fit_peaks_runs.py RUN_A RUN_B [--peaks]` | Cost-component comparison (hand-fit corpus: extents, sharing, continuum family). |
| `br_pinned.py RUN` | Lines the solve switched off with its branching-ratio nuisance, and the truth lines they cost. |
| `gate_rej.py RUN` | The keep gate's rejections in each problem's last plan, by predicted z: truth lines vs phantoms. |
| `refine_margins.py RUN...` | Refinement decisions: how close the score comparisons were, and rejected challengers that were healthier than the incumbent. |
| `worse_delivered.py RUN` | Final-solve ROIs that fit worse than no peaks but were delivered anyway, on the solve's own strongest-peak significance. |
| `fwhm_balloon.py RUN [RUN_B]` | How far each problem's accepted solve widths ran from the planner's width model, per energy band; with `RUN_B`, the problems whose low-energy ratio moved. |
| `side_extent.py HAND_RUN RUN [N]` | ROI side extents beyond the outermost peak, in truth FWHM, against the R500 hand fits (`run_hand.sh`). |
| `turnon.py RUN`, `turnon_rois.py RUN [v]` | Fitted peaks and delivered ROIs at the detector turn-on, classed by whether truth has a real line there. |
| `solve_health.py`, `roi_defects.py`, `fwhm_ratio.py`, `changed.py`, `hand_cmp.py`, `hand_split_rule.py` | Narrower diagnostics from earlier rounds. Each script's docstring says what it does. |

## Looking at fits
| Tool | What it does |
|---|---|
| `render_all.py RUN OUTDIR [--cmp OLD_RUN] [--ref]` | Whole-spectrum PNG per problem. New fit in red, old in blue, truth lines as green dash-dot lines, missed strong lines as orange triangles; `--ref` overlays the hand fit. |
| `overview.py RUN ID [ID...] --out DIR [--cmp OLD_RUN]` | The same images for chosen problems. |
| `plot_fit_peaks_rois.py RUN ID [--range lo,hi] [--worst N]` | Per-ROI zoom panels (data, background, fitted and reference continuum and total, truth lines, search peaks). |
| `triage_sum.py FILE...` | Verdict and defect-class tallies of a review (pass the part TSVs, e.g. `triage/c5/c5_part*.tsv`). |
| `review_page.py 'PART_GLOB' IMG_DIR OUT.html TITLE` | An HTML page of a review: spectra grouped by verdict, each with its notes and image. |
| `rr_compare.py 'OLD_GLOB' 'NEW_GLOB'` | Verdict transitions between a review and a later re-review. |
| `triage_page.py` | An older single-fit review page. |

## Fit-objective comparison (`peak_fit_objective_eval`)

Measures the statistics of `PeakFitLM`'s fit objectives - bias, pulls `(fit - reference)/reported
uncertainty`, coverage, failure and bad-fit-tail rates, CPU - running every objective on the identical
Poisson draw from identical starting values (paired), straight through `PeakFitLM::fit_peaks_in_roi_LM`
(no significance culls, so nothing is censored).  `ninja peak_fit_objective_eval`; `--help` for options.

- `--mode=synthetic`: spectra from InterSpec's own peak + continuum model (exact truth, no model
  mismatch): a primary grid of detector class x FWHM-in-channels x area x continuum level, plus
  continuum-type, skew and doublet suites (`--suites=`).
- `--mode=inject`: `peak_fit_accuracy_inject_compact` spectra; PCF record 3 is the noise-free
  spectrum, so `--replicas` draws are made from it.  Bias is reported against the GADRAS truth and
  against each objective's own pseudo-truth (its fit to the noise-free spectrum x1000), which removes
  peak-shape model mismatch from the statistical comparison; `pseudo_vs_chi2` shows how differently
  the objectives respond to that mismatch.  Low-resolution ROIs default to +-2.5 FWHM
  (`--lowres-roi-fwhm`), since wider ones take in structure a polynomial continuum cannot follow.
- Objectives (`--objectives=`): `chi2` (plain modified-Neyman, `NoSparseDataLikelihood`), `default`
  (no options: chi2 with sparse ROIs refit by Poisson likelihood), `likelihood` (every ROI refit,
  `ForcePoissonLikelihood`), and `chi2-cond` (chi2 with the old conditional area uncertainties - what
  the peak search and `fit_peaks_for_nuclides` use).  The likelihood fit is IRLS; to evaluate the
  all-in-Ceres alternative, build with `SPARSE_DATA_LIKELIHOOD_USE_CERES` set to 1 in
  `src/PeakFitLM.cpp`.  (The profiled-MLE and Mighell objectives of the 2026-10 evaluation were
  removed after it.)
- Outputs: `per_fit.tsv`, `summary.tsv` (per problem x peak x objective), `paired.tsv` (replicas where
  an objective differs from chi2 by > 3 sd), `run_meta.txt`; `--dump=SUBSTR:K` writes channel data and
  each objective's model for `tools/plot_objective_dump.py`.  `per_fit.tsv` also carries candidate
  sparse-data statistics of each fit (`sp_*`; `sp_peak` is the one PeakFitLM uses) and the number of
  ROIs the default fit refit by likelihood (`sparse_rois`).
- `tools/objective_report.py RUN --by area,cont_level --filter det=HPGe` tabulates (median over
  problems); `--paired` counts the tail; `--compare RUN_B` puts two runs side by side.
- Bit-identity of the default fit across a code change: diff the `chi2` rows of `per_fit.tsv`
  excluding the `cpu_s` column (`cut -f1-19,21`).
