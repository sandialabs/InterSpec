# Prompt: defect-driven improvement of `fit_peaks_for_nuclides`

How to use: fill in the **Scope** section, then paste everything below the line into a new session.  The
method was worked out improving scintillator (NaI) results in September 2026.  This is the generic
template; `HPGE_PROMPT.md` (the same effort for HPGe) and `LOWRES_FOLLOWON_PROMPT.md` (continuing the
scintillator work) in this directory are ready-made, tailored versions.  The tools are in this directory
(see `README.md`); the state that work reached is in the memory notes named at the end.

---

## Scope (edit before use)

- **Detector class / data**: e.g. *HPGe: the Detective-X inject corpus and the HPGe hand-fit corpus*, or
  *NaI: IdentiFINDER-R500-NaI, IdentiFINDER-NGH, SAM-Eagle-NaI-3x3 inject corpora*, or another detector
  directory under `$FPR_INJ` (LaBr3: Radseeker-LaBr3 and others; CZT: Kromek-GR1-CZT and others).
- **Must not move**: the other detector classes' output, bit for bit, unless we decide otherwise.
- **Branch / commits**: nothing is committed unless I ask; never list Claude as a co-author.

## Goal, and how to judge it

Improve `FitPeaksForNuclides::fit_peaks_for_nuclides` for the scope above.  Many people will use it every
day, and automated AI tools build on it, so it deserves sustained effort.  A good result is what an
experienced spectroscopist would accept: every visible line of the requested sources fitted, no peaks
where the data show none, ROIs drawn around features with sane continua, peak areas that match the data.

Work **defect-driven**.  Find *classes* of bad ROIs and fits, trace each class to the mechanism in the
code that produces it, fix the mechanism, and **look at the before/after fits**.  Aggregate truth scores
are guard rails, not the objective.  A change that raises matched counts but makes the fits look worse
is not progress.  A change that removes a real class of defect is progress even when the scores barely
move.

## Measurement kit

Everything is in `target/peak_fit_improve_ai/` (`README.md` there is the index): the harness in
`harness/`, the scripts in `tools/`, the reviewer rubrics in `review/`.  Scripts take their paths from
`tools/env.sh`: `FPR_WORK` (work area, default `~/fit_peaks_work`), `FPR_INJ` (inject corpus), `FPR_SITE` /
`FPR_DWELL` (inject site and dwell), `FPR_MANUAL` (hand fits), `FPR_DATA`, `FPR_PY` (a python3 with numpy
and matplotlib - the system python3 has no matplotlib: `python3 -m venv ~/fit_peaks_work/venv &&
~/fit_peaks_work/venv/bin/pip install numpy matplotlib`, then `export FPR_PY=~/fit_peaks_work/venv/bin/python`).

- **Build:** `cd target/testing/build_ninja && ninja fit_peaks_corpus_eval test_fitPeaksForSources test_RelActCalcAuto_ProfileApi`
  (the harness sources are in `target/peak_fit_improve_ai/harness/`, added to that build tree).
- **Snapshot:** run `snapshot.sh cNN` after every build you will measure, then measure only `$FPR_WORK/bins/eval_cNN`.
  Other sessions rebuild the shared build tree mid-run.
- **Corpora:**
  - `$FPR_INJ/<detector>/Livermore/300_seconds`: GADRAS-injected spectra with truth.
  - `$FPR_MANUAL`: hand fits. The HPGe set is the harness's default corpus; `run_hand.sh` runs 17 R500 NaI hand fits.
- **Harness flags that matter:**
  - `--det-type`: `Low`, `LaBr`, `CZT` or `High`. **The default is High**, so a scintillator run without it silently evaluates the HPGe configuration.
  - Always pass `--background=file --no-structure --weight=min_scored_energy=20`.
  - `--fit-timeout=1200` for scintillators, where U and Pu spectra take 2–5 minutes each.
  - `--problems a,b,c`, `--set field=value` (any registered config field; use it for sweeps), and `--debug ID` (one problem, serially, with the fitter's own trace).
- **Outputs of each run:**
  - `per_peak.tsv`: every truth and fitted peak with a verdict: matched, missed, extra, ghost, and so on.
  - `per_problem.tsv`: includes the solver's warnings, such as parameters pinned at a bound.
  - `roi_plan_trace.txt`: the planner's reasons. Every line group admitted or rejected, and why; every share or separate decision with its numbers; each ROI's continuum choice; a one-line summary of the first solve and of each refinement challenger.
  - `plot_data/*.json`: the data behind the images.
- **Tools:**
  - `cmp_runs.py OLD NEW`: strong/moderate/weak truth lines found, with the list of lines gained and lost. This is the main A/B tool.
  - `run_table.py PREFIX...`: found lines and extras per detector for `full.sh` runs, side by side.
  - `extras_cmp.py OLD NEW`: the extra peaks a candidate adds - read it whenever a change finds more lines.
  - `worse_delivered.py RUN`: ROIs delivered although they fit worse than no peaks.
  - `miss_why.py RUN minz maxz`: whether each miss was planned then lost, or rejected, and with what reason.
  - `invariants.py RUNs`: peaks outside their ROI, ghost peaks, timeouts.
  - `bitident.py`: per-peak identity between two runs.
  - `br_pinned.py RUN`: lines the solve switched off.
  - `render_all.py RUN OUT --cmp OLD`: images, new fit in red over old in blue.
  - `quick.sh`: a subset run plus images.
  - `full.sh BIN TAG [--set ...]`: all detectors plus the HPGe sets, in two parallel streams; `FPR_GUARD_REF` / `FPR_GUARD_SETS` for bit-identity.
  - `full_hpge.sh`: the HPGe sets (Detective-X 30/300/1800 s, Detective-EX, planar, Falcon, Fulcrum, hand fits).
  - `guard.sh` / `guard_lowres.sh`: the HPGe or the scintillator guard alone, with bit-identity.
  - `gallery.sh`: the harness's HTML galleries (`manual` = the HPGe hand fits).

## The loop

1. **Baseline.** Snapshot, run `full.sh`, build a gallery, and do a full visual review (below).  Record the
   numbers: strong/moderate/weak found, extras, raw cost per detector, and the review verdict counts.
2. **Triage.** Cluster the review's defect notes and the `miss_why.py` categories into classes.  Pick the
   class with the most, or the most severe, instances.
3. **Trace.** Take 2–3 exemplars.  Read their `roi_plan_trace.txt` sections, then run
   `--debug <id>` for the solve and refit detail, until you can name the exact code path.  A temporary
   env-var hook is fine while investigating, for example dumping the solver's initial parameters and
   bounds, which exposed two solver bugs.  **Remove it before any full run**: grep the diff for
   `getenv` and `TEMPORARY`.
4. **Fix the mechanism, gated.**
   - Add a `PeakFitForNuclideConfig` field that is off by default and on in the detector class's defaults
     (`s_default_non_hpge_config` for scintillators, which **also covers LaBr3 and CZT**).
   - Register it with `FPN_*_FIELD` so `--set` can switch it.
   - Changes to shared code (`PeakFitLM`, `RelActCalcAuto`) must be opt-in options, off by default.
     InterSpec's interactive fitting and the Isotopics tool use that code too.
   - Write the *why* in the comment and name the exemplar spectrum.
5. **Quick check.** Run `quick.sh` on the exemplars plus a regression set: every spectrum a previous
   review marked WORSE/MIXED/BAD, plus the exemplars of earlier fixes.  Run `cmp_runs.py` against the
   previous candidate, then **look at the images**, including those of spectra you did not target that
   changed.
6. **Full check.**
   - `full.sh` on every detector of the class.
   - `invariants.py` on every run.
   - `bitident.py` on the runs that must not move.
   - Both unit suites.
7. **Re-review** the changed and flagged spectra with the same batches and rubric, so the numbers stay
   comparable.
8. **Record.** Update the memory note and `TODO.md`, then take the next class.

## Visual review with agents

- Render with `render_all.py RUN OUTDIR --cmp PREVIOUS_RUN` and split the ids into batches of about 25.
  Launch one general-purpose agent per batch, 8 in parallel for 205 spectra, with the prompt below.  A
  batch costs 250–300k tokens, so after the first full review, re-review only flagged or changed spectra.
- Tally with `triage_sum.py` and the CMP counts.  `review_page.py` makes an HTML page; `rr_compare.py`
  compares a re-review with the review it follows.
- Reviewers have blind spots.  They call real but faint lines phantoms: escape peaks on a steep turn-on,
  blended lines.  **Check a "phantom" or "junk ROI" claim against `per_peak.tsv`** (is the fitted peak
  `matched`?) before acting on it.  The physics notes in `RUBRIC.md` exist because reviewers without them
  mis-called real iodine escape peaks.
- `RUBRIC.md` is written for NaI.  For HPGe, rewrite its numbers first: ROI widths in FWHM, sideband
  expectations, x-rays that resolve, only pair-production escape peaks.

The agent prompt used (substitute `<N>`, `<NEW>`, `<OLD>`, the image directory and the paths):

```
You are reviewing automated gamma-spectrum peak fits (<detector type> detector) by LOOKING at images, as
an expert spectroscopist would.
Directory with instructions: <repo>/target/peak_fit_improve_ai/review
1. Read `RUBRIC.md` there first (how to read the images, defect classes, verdicts) and view its three
   calibration example images.
2. Then read `RUBRIC_BOTH.md` there: each image shows TWO fits - RED = the new fit (run <NEW>), BLUE =
   the old fit (run <OLD>). You give an absolute verdict on RED plus a RED-vs-BLUE comparison.
Your spectra: the ids listed one per line in `<batch file>`.
Images: `<image dir>/<id>.png` (view each one with the Read tool; look carefully at every panel).
Write your results, exactly in the TSV format RUBRIC_BOTH.md specifies (tab-separated, no header; defect
lines then exactly one CMP line per spectrum), to `<out dir>/<NEW>_part<N>.tsv`. Append as you go (so
partial progress is kept). Do every spectrum in your list. Judge only what you see in the images; do not
read or reason about the software source code. When finished, reply with a two-line summary: counts of
GOOD/MINOR/BAD/CATASTROPHIC and of BETTER/WORSE/MIXED/SAME.
```

## Rules learned the hard way

- **Noise level.** The rel-eff solve is non-convex, so small plan changes flip marginal outcomes.  A
  shift of ±2–4 strong lines per detector between candidates is noise, and so is a single timeout.
  Judge by mechanisms and exemplars.  Read the gained/lost lists, not the totals.  When a strong line
  flips, re-run it on both binaries with `--debug` and find the step that changed.
- **Determinism.** The same binary gives identical output on a subset run and a full run.  When two
  runs disagree, it is a real code difference (or a timeout), never randomness.
- **Run every detector after a planner gate change.**  A data test tuned on R500 (2.9 keV/channel)
  broke SAM (12.5 keV/channel, where ±2 FWHM is about 5 channels) and NGH (a neighbour line's flank, and
  a background Cs137 line inside the test window).
- **Invariants catch what scores hide.**  Fitted peaks outside their own ROI and zero-amplitude "ghost"
  peaks were found by `invariants.py`, not by the truth score.
- **Fixes that worked on their exemplar and broke something else** (all caught only by looking wider):
  - Gluing escape groups to neighbours made 5–13 FWHM ROIs.
  - Floating peaks for x-rays exposed a solver seeding bug.
  - Exempting escape peaks from a filter produced ghost peaks.
  - A per-side "room" rule glued a 15–132 keV chain together through a weak group.
- **Compare truth-line counts, not raw costs, across scorer changes.**  Keep the truth denominator fixed.
- **HPGe (or whatever must not move) guard:** per-peak `bitident.py`.  Equal summary lines can hide
  compensating differences.
- **Mention review prompt changes when reporting**, for example added physics notes: they change what
  reviewers call a phantom.
- **Truth quirks** (GADRAS):
  - The truth includes detector effects: iodine escape peaks (`EscapeXRay` components) and pair-production escapes.
  - Truth energies can be intensity-weighted centroids of blends.
  - Some truth lines disagree with the decay data (Np237's 405 keV line at 0.58 of 312 keV vs 0.08 in the data; I123/In111 K x-rays roughly 10⁴–10⁵ above the decay-data yields).
  - Before calling something a fitter error, check whether it is a data discrepancy.
- **Hygiene:** never delete user data files, kill only your own processes, and never broad-`pkill`.
- **A refactor must leave the default path bit-for-bit alone.**  Before a full run, re-run 3-5 problems
  and compare the exact numbers the trace prints (refinement scores, solve summaries) with the
  reference run.  A flag that already meant something (`delivered_model` chose the significance test
  AND was reused for a new selection) silently changed every NaI refinement for two full runs.
- **The refinement comparison is where changes turn into lotteries.**  Any change to the plans or the
  score reshuffles which pass is accepted, and strong lines move both ways.  Three scores were
  measured (Sept 2026): common-domain (default, and a corrected variant), whole-range against the
  SNIP continuum, with and without continuing past a rejected pass.  The whole-range score found +16
  strong / +33 moderate lines but a visual review called it worse twice as often as better: on NaI
  the SNIP continuum does not follow scatter / backscatter / Compton-edge structure, so an ROI
  modelling that structure with peaks scores as progress.  More truth lines is not the goal - look.
- **A CZT's discriminator slices lines**: the GR1 is dead below ~39 keV; a line centred there shows as
  a narrow half-peak, the search finds it, and it can set the whole width model.  Check where the
  first live channel is before trusting any low-energy peak.

## Mechanisms found so far (recognise them quickly)

- **Width model.** The sqrt-polynomial fit collapsed to a constant FWHM (a monotonic prefilter from
  0 keV).  It was replaced by a class-prior shape times a robust scale k (`detail::fit_fwhm_function_robust`;
  NaI k clamped to 0.8–1.7).  The planner's width now seeds the solver (`starting_fwhm_*`).
- **Refinement acceptance** once rewarded a broken ROI that the final filter later removed.  The
  challenger and the incumbent are now scored on their common segments, with anchor guards.
- **Solve collapses** (activity driven to ~0, or an ROI worse than no peaks):
  - Lines outside every ROI extrapolated at absurd efficiency → `solve_lines_within_roi_span`.
  - Decay-data x-ray yields orders of magnitude low → the zero-activity retry also fires when a
    nuclide's *principal anchor* (highest-yield matched peak) is explained < 10%.
  - The branching-ratio nuisance could zero a line for (1/σ)² ≈ 27 chi2 → `rel_eff_auto_br_min_yield_fraction`.
  - RelActCalcAuto's initial estimate matches search peaks to lines of any yield (open; see `TODO.md`).
- **Observable refit** (`compute_observable_peaks`):
  - PeakFitLM's shared width model keeps the last peak within ±20% of the first → `observable_independent_widths`.
  - The ROI edge moved whenever the edge peak's mean shifted → `observable_edge_moves_on_removal_only`.
  - A collapse guard reverted to a solve that was itself broken → it now requires the solve to have followed the data.
  - Escape companions are held at their tied amplitude.
- **ROI geometry on coarse binning.**  Leak-based splitting is sound, but it made 3-channel ROIs
  (`share_min_side_channels`, `share_min_roi_channels`), and the one-channel gap was always taken from
  the upper ROI (`roi_gap_at_boundary`).
- **Admission gates** (`plan_rois_impl`):
  - The data tests use `detail::fixed_shape_peak_z`, a Gaussian over a quadratic, with neighbour lines, net of the background.
  - SNIP is unreliable on the detector turn-on.
  - Near the analysis floor, windows must be clipped.
  - An escape group is judged by area consistency with its parents.
  - Background-coincident lines fail the net test (La138 1436 vs K40 1461; open).

## Code map

- `src/FitPeaksForNuclides.cpp`:
  - `plan_rois_impl`: the planner.
    - Line groups and admission (keep gate, refutation, data tests, data-detected/data-evident, sub-extent rules, escape groups).
    - Sharing (leak cut, valley, room rules), span cap, merges.
    - Extents (visible-line core, `extend_roi_by_sidebands`), clipping, floor at the extent.
    - Continuum type, channel alignment; every decision is written to the trace.
  - `cluster_gammas_to_rois`: turns predicted lines into planner input.
  - `estimate_initial_rois_using_relactmanual` / `_fallback`: the first plan, from a manual rel-eff fit on search peaks.
  - `fit_peaks_for_nuclide_relactauto`:
    - The first solve and its retries (no-ecal, physical-model desperation, zero-activity retry).
    - The refinement loop (re-plan, challenger vs incumbent, anchor guards), the final filter, combining peaks.
    - Observable peaks, labels and result assembly.
  - `compute_observable_peaks`: the per-ROI PeakFitLM refit on raw data (initial filter, held escapes, neighbour tails, refit loop and edge rule, edge-dive re-measure, collapse guard).
  - `detail::fit_fwhm_function_robust`, `detail::fixed_shape_peak_z`.
  - `PeakFitForNuclideConfig` is declared in `InterSpec/FitPeaksForNuclides.h`, with defaults in `s_default_non_hpge_config` and the HPGe default, registered through `FPN_*_FIELD`.
- `src/RelActCalcAuto.cpp`: the joint solve.
  - `Options` (`iodine_escape_peaks`, `model_lines_outside_roi_span`, `additional_br_uncert`, `additional_br_min_yield_fraction`, `starting_fwhm_*`).
  - `peaks_for_energy_range_imp` builds each ROI's model peaks, including escape companions.
  - The BR-nuisance parameters, the initial estimate, and the Ceres setup.
- `src/PeakFitLM.cpp`: the LM refit and its width parametrisations.
- `target/peak_fit_improve_ai/harness/`: `fit_peaks_corpus_eval.cpp`, `FitPeaksCorpusScore.*` and `FitPeaksCorpusReport.*` - the harness, its verdicts and its galleries.

## Lessons from the September 2026 scintillator rounds

- **Planner geometry changes re-decide the refinement passes.** Any change to where the planner puts
  ROIs (a keep gate, an extension cap, a core width) flips near-tie refinement decisions, and strong
  lines move both ways on unrelated spectra.  Fix presentation problems (ROI width, joins) DOWNSTREAM
  of the solve, in the observable refit, where they cannot feed back.
- **A rule that deletes a phantom ROI can break its neighbour.** The deleted ROI's boundary may have
  been what kept the next ROI out of a region the solve cannot model (Lu177m_Sh).
- **Measure a change alone before stacking it.** Several of the round's reverted changes looked right
  on the spectra they were written for and cost more elsewhere.
- **The truth score cannot see a broken fit** whose peaks happen to sit near the truth energies; the
  visual review found the ones the numbers called gains.

## HPGe-specific notes

- **Reference:** the hand-fit corpus scores ROI sharing and extents against human fits, as well as
  truth.  Use it together with the Detective-X inject truth.
- **The guard runs the other way:** the scintillator runs of the current candidate must stay
  bit-identical.
- **Resolution is ~10–20 channels per FWHM.** Per-channel rules from the scintillator work do not apply,
  and x-ray lines resolve.  Iodine escapes do not exist; pair-production SE/DE escapes above 1.6 MeV are
  handled by `add_escape_peak_floating_peaks_if_appropriate`.
- **State on 2026-09-06** is in memory note `project_fit_peaks_hpge_state.md`.

## Where things stand

- **The code is on `upgrade/Wt4`**, in the commit "Agentic peak fit experiment to improve low-res fits"
  (2026-09-25; the 2026-09-26 scintillator round was squashed into it on 09-27).  The local branches
  `tmp/fit-peaks-overhaul` and `review/fit-peaks-overhaul` keep the earlier step history.
- **Memory notes** (the Claude Code project memory for this repository, `~/.claude/projects/<repo-path>/memory/`):
  - `project_fit_peaks_nai_state.md`: every scintillator round, its numbers, and every mechanism above
    with exemplars.
  - `project_fit_peaks_hpge_state.md`.
  - `feedback_fit_peaks_defect_driven.md`.
  - `feedback_fit_peaks_generality.md`: stay general, with no DRF or prior information assumed.
- **`TODO.md`** at the repository root: the sections "fit_peaks_for_nuclides on scintillators",
  "fit_peaks_for_nuclides merged onto upgrade/Wt4 - follow-ups" and "fit_peaks_for_nuclides
  scintillator round 2026-09-26".
- **Scintillator numbers now** (run `c159`; strong / moderate truth lines found, truth lines >= 20 keV
  that are not don't-care): R500 96.7 / 52.6 %, NGH 95.4 / 43.9 %, SAM 89.8 / 44.5 %, LaBr3
  96.9 / 54.7 %, CZT 93.8 / 35.6 %.
- **Reference runs** in `~/fit_peaks_work/runs/` (see `~/fit_peaks_work/NOTES.txt`): `c159_<set>`, the
  code as committed now, for `<set>` in r500 ngh sam labr czt detx manual.  The matching binary is
  `~/fit_peaks_work/bins/eval_c159`; galleries are in `~/fit_peaks_work/galleries/c159/` (scintillators)
  and `~/fit_peaks_work/galleries/mergedB/` (HPGe: `c159_detx` / `c159_manual` are bit-identical to
  `mergedB_detx` / `mergedB_manual`).  The older scintillator runs `mergedB_<set>` (the 09-25
  landing), `c77_<set>` and `ref0924_<set>` are history.

## Reporting

After each candidate, report to me:

- A numbers table per detector (strong, moderate, weak, extras, raw cost), against the previous candidate and the baseline.
- The visual review counts.
- Each change, the mechanism it fixed, and its exemplar.
- What got worse, and why.
- The open items.

Keep going until the defect classes run out or I stop you.
