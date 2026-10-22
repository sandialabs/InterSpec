# Prompt: improve `fit_peaks_for_nuclides` on HPGe - defect-driven, judged by looking

How to use: paste everything below the line into a new Claude Code session in this repository.  It
repeats, for HPGe, the method used on the scintillators (NaI, LaBr3, CZT) in September 2026.  Update
the "Where HPGe stands" section if much time has passed.

---

## Goal

Improve `FitPeaksForNuclides::fit_peaks_for_nuclides` on HPGe spectra.  It is the default way users fit
peaks for sources in InterSpec and the workhorse behind the LLM assistant, so it deserves sustained,
careful effort.  A good result is what an experienced HPGe spectroscopist would accept: every visible
line of the requested sources fitted, no peaks where the data show none, ROIs drawn around features
with sane continua, peak areas that match the data.

Work **defect-driven and judge by looking.**  Find *classes* of bad ROIs and fits, trace each class to
the mechanism in the code that produces it, fix the mechanism, and look at the before/after fits.
Truth-matched counts are guard rails, not the objective: in the scintillator round a scoring change
found 51 more truth lines and a visual review still called it worse twice as often as better.  A change
that removes a real class of defect is progress even when the numbers barely move.

Standing rules from the user:
- Nothing is committed unless I ask, and never list Claude as a co-author.
- Project TODOs go in `TODO.md` at the repository root (git-ignored).
- Never delete user data files.  Kill only processes you started; never broad-`pkill`.
- Follow `CLAUDE.md`'s code style.
- The fitter must stay general: no DRF or other prior information assumed, low to high statistics, any
  detector.  Prefer physically motivated, dimensionless criteria over thresholds tuned to one corpus.

## Where HPGe stands

- **HPGe was the first target** (2026-09-05 to 09-07).  That work built the single-pass ROI planner,
  the physics envelope (sibling-absence checks), data-detected admission, escape peaks taken from the
  data, the step-continuum rules and the ROI core and extent rules.  On 2026-09-06 it found these
  shares of strong truth lines (fitter only):

  | Detector | Strong lines found |
  |---|---|
  | Detective-X | 99.5 % |
  | Detective-EX | 98.7 % |
  | HPGe planar 50 % | 97.6 % |
  | Fulcrum 40h | 89.9 % |

  The 168 hand-fit spectra sat at 92.3 % of reference peaks above z=3.  Fulcrum 40h has only 37 usable
  problems: 186 of its 223 truth files are empty.
- **HPGe fits have never had a visual review.**  The review-by-agents loop was developed later, on
  NaI.  Expect ROI, continuum and presentation defects the truth score cannot see.  The scintillator
  reviews found, among others:
  - over-wide ROIs, and ROIs joining features across a valley;
  - continua diving to zero at an ROI edge, with an edge peak compensating;
  - comb peaks on smooth humps, and phantom peaks on edges;
  - ROIs on the detector turn-on.
- **HPGe output was held bit-identical through the scintillator rounds** (2026-09-22 to 09-24, and
  2026-09-26): every scintillator change was gated off HPGe.
- **The code is on `upgrade/Wt4`**, in the single commit "Agentic peak fit experiment to improve
  low-res fits" (2026-09-25; the 2026-09-26 scintillator round was squashed into it on 09-27).  The
  local branches `tmp/fit-peaks-overhaul` and `review/fit-peaks-overhaul` keep the step history.
  The 09-25 landing moved HPGe slightly, through
  upgrade/Wt4's peak-CDF step rework (PeakFitLM now fits the step coefficient) and solver changes,
  and through the NoisePlusCurvedPower FWHM form, which every class's solves now use (review-branch
  commit 2fa84c17; `git show 2fa84c17 | git apply -R` removes it while that branch exists).
  Strong / moderate / weak truth lines found: Detective-X inject 991/555/71 at `ref0924`, 991/559/69
  as landed; the hand fits 920/470/14 and 920/469/14.
- **References**: `~/fit_peaks_work/runs/mergedB_detx` and `mergedB_manual` describe the HPGe state as
  landed (binary `bins/eval_mergedB0`, galleries in `~/fit_peaks_work/galleries/mergedB/`);
  `c159_detx` / `c159_manual`, from the code as committed now, are bit-identical to them.
  `c159_<set>` is the scintillator guard (binary `bins/eval_c159`).  `ref0924_*` is history.
  `~/fit_peaks_work/NOTES.txt` says what each saved item is, and `TODO.md`'s section
  "fit_peaks_for_nuclides merged onto upgrade/Wt4 - follow-ups" lists what the merge surfaced.
- **Known open HPGe items** (memory note `project_fit_peaks_hpge_state.md`, plus `TODO.md`):
  - The step-continuum gate (`step_cont_min_peak_significance` 40) was tuned on Detective-X and misses
    visible steps at z 21-27 on other detectors.
  - ROIs are sometimes clipped to ~0.7 FWHM on one side after planning, which also makes continua look
    unphysical.
  - The automated search finds too few peaks on some spectra (8 in a whole Fulcrum U235 spectrum).
  - RelActCalcAuto's internal 4-FWHM gain freedom disagrees with the outer 2-FWHM drift bound.
  - The Pb-for-shielded rule is over-broad.
  - RelActCalcAuto sometimes returns its own seed; U235_Unsh_5000 loses its z=50 185.7 keV line this
    way.
  - The R6 interferer logic can hide a genuine source line (Ag110m 687 keV taken as Ra226).
  - Self-fluorescence x-rays are rejected by the envelope (Lu176).
- **Scintillator-round findings that apply to HPGe and were never tried there.**  Measure each as an
  early, cheap experiment; they are all `--set`-able:
  - `PeakFit::chi2_for_region` ignores the channel range it is given whenever the continuum has an
    energy range (its locals shadow the parameters).  The HPGe refinement score
    (`solution_chi2_over_segments`) relies on it, so each ROI is charged its whole-range chi2 once per
    shared segment.  Taking the peaks *centred* in an ROI also pulls in a neighbour's copy of a line
    along with the neighbour's continuum.  This biases against splitting joined ROIs.  The corrected
    score is `refinement_delivered_segments=1`, which needs `roi_significance_delivered_model=1`.
    That second flag, the ROI significance test on the delivered model, is itself off for HPGe
    "until measured there".
  - **The whole-range refinement score** (`refinement_whole_range_score=1`, best with
    `refinement_keep_best=1`, also needs `roi_significance_delivered_model=1`) scores every candidate
    over the whole analysis range: significant ROIs by their own model, everything else by the shared
    SNIP continuum.
    - It fixes the "refinement lottery": the loop otherwise stops at the first rejected pass, and the
      common-domain score cannot credit the lines a challenger adds.
    - It failed on NaI only because the SNIP continuum does not follow NaI's broad scatter,
      backscatter and Compton-edge structure.
    - On HPGe, with 1-2 keV peaks on a continuum that is smooth at that scale, the SNIP should follow
      it well.  **This is the most promising first experiment.**
  - `tools/worse_delivered.py RUN` counts ROIs that fit worse than no peaks yet were delivered on the
    solve's own strongest-peak significance.  On scintillators that was 11-35 per detector; count it
    on HPGe.
  - The manual rel-eff ladder forgives lines a form cannot predict, which lost NGH Ca47_Sh's z=53
    1297 keV line.  The code is shared with HPGe.

## Measurement kit

Everything is in `target/peak_fit_improve_ai/` (`README.md` there is the index): the harness in
`harness/`, the scripts in `tools/` (they source `tools/env.sh`), the reviewer rubrics in `review/`.

- **Setup, once.**
  - `python3 -m venv ~/fit_peaks_work/venv && ~/fit_peaks_work/venv/bin/pip install numpy matplotlib`,
    then `export FPR_PY=~/fit_peaks_work/venv/bin/python`.  The system python3 has no matplotlib.
  - The work area defaults to `~/fit_peaks_work`.
- **Build:** `cd target/testing/build_ninja && ninja fit_peaks_corpus_eval test_fitPeaksForSources test_RelActCalcAuto_ProfileApi test_fitPeaksCorpusScore`.
- **Snapshot** every build you will measure: `tools/snapshot.sh hN`, then measure only
  `~/fit_peaks_work/bins/eval_hN`.  Other Claude sessions share `build_ninja` and rebuild it mid-run.
  Use tags `h0`, `h1`, ... so HPGe runs never collide with the scintillator runs (`c159`, `mergedB`, ...).
- **HPGe corpora:**
  - **The hand fits**: 168 Detective-X 300 s spectra, the harness's default corpus.
    - `tools/env.sh`'s `fpr_hand_run`, or `$BIN --datadir=$FPR_DATA --out=DIR --background=file`.
    - They are scored against the hand fits for presence AND structure (ROI sharing, extent sides,
      continuum family), with the GADRAS-inject truth attached.
    - Peaks the hand fits omitted are often real, so settle ambiguity with the truth.
  - **GADRAS-inject** (`--corpus-format=inject --no-structure`, `--det-type=High`, which is the default):
    - Detective-X at 30 / 300 / 1800 s (statistics robustness);
    - Detective-EX, HPGe_Planar_50%, Falcon 5000 (not yet used), Fulcrum40h;
    - Detective-X_noskew, by its name the same detector simulated without peak skew (useful to
      isolate skew handling; confirm before relying on it).
    - Baltimore and Denver backgrounds exist besides Livermore (`FPR_SITE`).
  - `tools/full_hpge.sh BIN TAG [--set ...]` runs all of these in two parallel streams;
    `tools/quick.sh` runs a subset with images; `--problems a,b,c` and `--debug ID` (one problem,
    serially, with the fitter's own trace) work on any run.
- **Outputs of each run:**
  - `per_peak.tsv`: every truth, reference and fitted peak with its verdict.
  - `per_problem.tsv` and `per_roi.tsv`.
  - `roi_plan_trace.txt`: the planner's reasons for every line group, share or separate decision and
    continuum choice, plus a one-line summary of the first solve and of each refinement challenger with
    its score.
  - `plot_data/*.json`, the data behind the images.
- **Scorer verdicts to know:**
  - `annih`: a fitted 511 keV peak, free.
  - `dontcare_nodata`: a truth line the data do not show.
  - `dontcare_shield_xray`: Pb K x-rays no requested source explains.
  - `legit`: an unmatched fitted peak on a real truth photopeak, free.
  - `extra_bkg`: a fitted peak on a background line.
  - Compare `matched` counts on a fixed denominator, never raw costs across scorer changes.
- **Compare with:**
  - `cmp_runs.py OLD NEW`: strong/moderate/weak truth lines found, plus the gained/lost list.
  - `compare_fit_peaks_runs.py OLD NEW --peaks`: the hand-fit cost components.
  - `extras_cmp.py OLD NEW`, `miss_why.py`, `invariants.py`, `worse_delivered.py`, `refine_margins.py`.
  - `changed.py OLD NEW`: the problems whose fitted ROIs or peaks changed.  The lines starting with a
    letter are the ids; the last line is the count.
- **Guard: the scintillators must stay bit-identical.**  Every scintillator fit uses
  `s_default_non_hpge_config`; HPGe uses the plain header defaults (`PeakFitForNuclideConfig::default_config`).
  - For an HPGe-only switch, default the new field ON in `InterSpec/FitPeaksForNuclides.h` and set it
    OFF in `s_default_non_hpge_config`, or gate the code on `det_type == High`.
  - Register every new field with `FPN_*_FIELD` so `--set` can flip it.
  - Changes to shared code (`PeakFitLM`, `RelActCalcAuto`, `PeakFit`) must be opt-in options, off by
    default: InterSpec's interactive fitting and the Isotopics tool use that code too.
  - At the end of each class run `tools/guard_lowres.sh BIN TAG c159` (about 1.5 h).  It checks
    all five scintillator sets peak by peak against the saved reference runs.  If upgrade/Wt4 has
    moved the scintillator results since the landing, make your own low-res baseline and guard
    against that instead.

## Visual review for HPGe: build this first

The image and rubric tools were made for NaI, where a whole spectrum fits in three panels.  An HPGe
spectrum has 5-40 ROIs, each a few keV wide, so first build the HPGe versions:

1. **An HPGe review image.**  One PNG per spectrum:
   - a small log-scale whole-spectrum strip with the ROI bands;
   - a grid of per-ROI zoom panels, each the ROI ± ~2 of its widths, linear y;
   - in each panel: data, the new fit in red (continuum dashed, total solid, peak marks with z), the
     old fit in blue (`--cmp`), the hand fit in green (`--ref`, for the 168 hand-fit spectra), truth
     lines, and missed strong truth lines;
   - panels for strong truth lines that got no ROI, so misses are visible too.

   `tools/plot_fit_peaks_rois.py` already draws per-ROI panels with the reference and truth; extend it
   with `--cmp`, or add a grid mode to `tools/overview.py`.  Keep it in `tools/`.
2. **`review/RUBRIC_HPGE.md`**, derived from `review/RUBRIC.md`.  Keep its defect codes, verdicts and
   TSV format, so `triage_sum.py`, `review_page.py` and `rr_compare.py` keep working.
   - **Rewrite the expectations for HPGe.**  Measure what the hand fits do: their typical ROI widths
     and side extents in FWHM, and when they share an ROI between neighbours.
   - **Physics notes the reviewers need**, so they do not mis-call real peaks:
     - x-ray doublets resolve (K-alpha1/alpha2, K-beta);
     - step continua under strong peaks;
     - low-energy tailing (skew) of HPGe peaks;
     - pair-production single and double escapes, 511/1022 keV below lines above 1.6 MeV;
     - Ge K x-ray escape peaks ~10-11 keV below strong low-energy lines;
     - sum peaks;
     - Pb K x-rays from shielding;
     - a broad backscatter hump near 200-250 keV and Compton edges, which are not peaks;
     - the Doppler-broadened 511 keV annihilation peak, wider than its neighbours;
     - truth markers that can be blend centroids.
   - **Calibration.**  Pick three spectra (GOOD, MINOR/BAD, CATASTROPHIC or BAD), render them into
     `review/calibration_hpge/`, and write their expected verdicts and defects into the rubric.  This
     is what makes different agents grade alike.
3. **Review with agents.**
   - One general-purpose agent per batch of ~20 spectra, 5-8 in parallel, each told to append its
     TSV as it goes; a batch costs 250-300k tokens.
   - First a baseline review of all 168 hand-fit spectra, plus ~50 from Detective-EX / planar /
     Falcon / Fulcrum.
   - Later, re-review only changed (`changed.py`) and previously flagged spectra, with the same
     batches and rubric, so verdict counts stay comparable.
   - Tally with `triage_sum.py FILE...`; `review_page.py` makes an HTML page to hand to me.
   - Spot-check reviewer claims yourself by viewing the images.  Before acting on a "phantom" or
     "junk ROI" call, check it against `per_peak.tsv` (is the peak `matched` or `legit`?): reviewers
     call real but faint lines phantoms.

   The agent prompt (substitute the paths, `<N>`, `<NEW>`, `<OLD>`):

   ```
   You are reviewing automated gamma-spectrum peak fits (HPGe detector) by LOOKING at images, as an
   expert spectroscopist would.  Instructions: <repo>/target/peak_fit_improve_ai/review.
   1. Read `RUBRIC_HPGE.md` there first and view its calibration images.
   2. Then read `RUBRIC_BOTH.md`: each image shows the new fit (RED, run <NEW>) and the old fit (BLUE,
      run <OLD>); give an absolute verdict on RED plus a RED-vs-BLUE comparison.
   Your spectra: the ids listed one per line in `<batch file>`.  Images: `<image dir>/<id>.png` (view
   each with the Read tool; look at every panel).  Append your results, in the TSV format the rubric
   specifies, to `<out dir>/<NEW>_part<N>.tsv` as you go.  Do every spectrum.  Judge only what you see;
   do not read the software source.  When finished, reply with two lines: counts of
   GOOD/MINOR/BAD/CATASTROPHIC and of BETTER/WORSE/MIXED/SAME.
   ```

## The loop

1. **Baseline.**
   - Snapshot `h0` and check that it is bit-identical to `mergedB_detx` / `mergedB_manual` (equal to
     `c159_detx` / `c159_manual`).  If it is not, something on upgrade/Wt4 has changed the fitter
     since 2026-09-27: note what, and use `h0` as the reference.
   - Run `full_hpge.sh`, build the galleries (`gallery.sh BIN h0 Detective-X:High manual`), and do the
     baseline visual review.
   - Record the numbers per set (strong/moderate/weak found and extras for inject; raw cost and
     components for the hand fits) and the review verdict counts.
2. **Cheap experiments first**: the flag combinations above, each via `--set` on `full_hpge.sh`,
   each with a visual review of the changed spectra.
3. **Triage.**  Cluster the review's defect notes and `miss_why.py` reasons into classes.  Take the
   class with the most, or the most severe, instances.
4. **Trace** 2-3 exemplars through `roi_plan_trace.txt`, then `--debug <id>`, until you can name the
   exact code path.  Temporary env-var hooks are fine while investigating; remove them before any
   full run (grep the diff for `getenv` and `TEMPORARY`).
5. **Fix the mechanism, gated** (see the guard), with the *why* and the exemplar spectrum written in
   the code comment.
6. **Quick check** the exemplars plus a regression set (every spectrum a review marked WORSE, MIXED
   or BAD, and earlier fixes' exemplars), with `quick.sh` and `FPR_CMP`.  Look at the images.
7. **Full check**:
   - `full_hpge.sh` and `invariants.py`;
   - `guard_lowres.sh` against `c159` (or your low-res baseline);
   - the unit suites (`test_fitPeaksForSources`, `test_RelActCalcAuto_ProfileApi`,
     `test_fitPeaksCorpusScore`; add a unit test for a new mechanism);
   - re-review the changed spectra.
8. **Record** in the memory note `project_fit_peaks_hpge_state.md` and in `TODO.md`, then take the
   next class.

## Rules learned the hard way (September 2026, scintillators; they carry over)

- **Noise.**  The rel-eff solve is non-convex, and small plan changes flip which refinement pass is
  accepted.  A shift of ±2-4 strong lines per set between candidates is noise.  Read the gained/lost
  lists, not the totals.  When a strong line flips, re-run it on both binaries with `--debug` and find
  the step that changed.
- **Determinism.**  The same binary gives identical output on a subset and on a full run.  A
  difference is always a code difference (or a timeout).
- **A refactor must leave the default path bit-for-bit alone.**  Before a long run, re-run 3-5
  problems and compare the exact numbers the trace prints (refinement scores, solve summaries) with
  the reference run.  In the last round, a flag that already chose one thing was reused for a second
  thing, and every NaI refinement silently changed for two full runs.
- **Measure a change alone before stacking it**, and run every set after a planner gate change: a rule
  tuned on one detector broke others.
- **Planner geometry changes re-decide the refinement passes.**  Fix presentation (ROI width, joins)
  downstream of the solve, in the observable refit, where it cannot feed back.
- **A rule that deletes a phantom ROI can break its neighbour**, whose boundary it had been holding.
- **Invariants catch what scores hide** (`invariants.py`: peaks outside their own ROI, ghost peaks,
  timeouts).
- **Truth quirks** (GADRAS):
  - Truth includes pair-production escapes.
  - Energies can be intensity-weighted blend centroids.
  - Some truth lines disagree with the decay data.
  - Source mixtures follow the mapping in `load_inject_problem` (U233 → U, U232, U233, and so on).
  - Before calling something a fitter error, check whether it is a data discrepancy.
- **InterSpec labels lines by the ultimate parent** (U238 owns 1001 keV).  Never "correct" that.

## Code map

- `src/FitPeaksForNuclides.cpp`:
  - `plan_rois_impl` (the planner): line groups, admission, sharing, span cap, extents, clipping,
    continuum type.  Every decision is written to the trace.
  - `estimate_initial_rois_using_relactmanual`: the first plan, from a manual rel-eff fit.
  - `fit_peaks_for_nuclide_relactauto`: the first solve and its retries, the refinement loop
    (`score_pair`, the anchor guards), the final filter, observable peaks and labels.
  - `compute_observable_peaks`: the per-ROI PeakFitLM refit on raw data.
  - `detail::fit_fwhm_function_robust`: the width model.
  - `add_escape_peak_floating_peaks_if_appropriate`: escape peaks.
- `InterSpec/FitPeaksForNuclides.h`: `PeakFitForNuclideConfig`, with every field documented;
  `default_config()` is in the .cpp.
- `src/RelActCalcAuto.cpp`: the joint solve.  `src/PeakFitLM.cpp`: the LM refit.
  `src/PeakFit.cpp`: `chi2_for_region`, the search.
- `target/peak_fit_improve_ai/harness/`: the harness, its verdicts and galleries.

## Where to read first

- `target/peak_fit_improve_ai/README.md`.
- Memory notes:
  - `project_fit_peaks_hpge_state.md`, `project_fit_peaks_corpus_harness.md`,
    `project_fit_peaks_escape_peaks.md`.
  - `feedback_fit_peaks_defect_driven.md`, `feedback_fit_peaks_generality.md`.
  - `project_fit_peaks_nai_state.md`, the scintillator rounds and their lessons.
- `TODO.md`, section "fit_peaks_for_nuclides on scintillators", for the shared items.

## Reporting

After each candidate, report to me:
- a numbers table per set against the previous candidate and the baseline;
- the visual review counts;
- each change, the mechanism it fixed, and its exemplar;
- what got worse, and why;
- the open items.

When a gallery or review page is ready, `open` it for me.  Keep going until the defect classes run
out or I stop you.
