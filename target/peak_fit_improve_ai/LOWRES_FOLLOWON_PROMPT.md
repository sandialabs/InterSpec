# Prompt: continue improving `fit_peaks_for_nuclides` on scintillators (NaI, LaBr3, CZT)

How to use: paste everything below the line into a new Claude Code session in this repository.  It
picks up where the September 2026 scintillator rounds stopped.  If more work has happened since,
update "Where things stand" first (`git log`, the memory note `project_fit_peaks_nai_state.md`).

---

## Goal

Keep improving `FitPeaksForNuclides::fit_peaks_for_nuclides` on low-resolution detectors: NaI/CsI,
LaBr3, CZT.  It is the default way users fit peaks for sources in InterSpec and the workhorse behind
the LLM assistant, so it deserves sustained, careful effort.  A good result is what an experienced
spectroscopist would accept: every visible line of the requested sources fitted, no peaks where the
data show none, ROIs drawn around features with sane continua, peak areas that match the data.

Work **defect-driven and judge by looking.**  Find *classes* of bad ROIs and fits, trace each to the
mechanism that produces it, fix the mechanism, and look at the before/after fits.  Truth-matched
counts are guard rails, not the objective.  The last round showed why: a refinement score that found
51 more truth lines was called worse about twice as often as better by the visual review.

Standing rules from the user:
- Nothing is committed unless I ask, and never list Claude as a co-author.
- Project TODOs go in `TODO.md` at the repository root (git-ignored).
- Never delete user data files.  Kill only processes you started; never broad-`pkill`.
- Follow `CLAUDE.md`'s code style.
- The fitter must stay general: no DRF or prior information assumed, low to high statistics, any
  detector.
- **HPGe output must stay bit-identical** unless we decide otherwise.

## Where things stand (2026-09-27: the 2026-09-26 round landed on upgrade/Wt4)

- **The code is on `upgrade/Wt4`**, in the single commit "Agentic peak fit experiment to improve
  low-res fits".  It first landed on 2026-09-25; the 2026-09-26 round was squashed into it on
  2026-09-27 (the local branch `tmp/fit-peaks-lowres3` points at the same commit).  That round's steps
  are not commits: `TODO.md`'s section "fit_peaks_for_nuclides scintillator round 2026-09-26" and the
  memory note record each one.  The local branches `tmp/fit-peaks-overhaul` and
  `review/fit-peaks-overhaul` keep the earlier step history that the hashes in this prompt refer to:
  - on `tmp/fit-peaks-overhaul`: the HPGe work ends at 1a328a8d (2026-09-06); the scintillator
    rounds are d7b9ef81, 0d4857fd and 2e9db453 (run tag `c77`), then 74dea9e4 (the 2026-09-24
    round) and d0feec70 (tools moved);
  - on `review/fit-peaks-overhaul`: d1608001 (all of the above merged onto upgrade/Wt4), 23689c4b
    (RelActCalcAuto holds bound-pinned mass-fraction totals and absent nuclides out of the rank
    analysis) and 2fa84c17 (the NoisePlusCurvedPower FWHM form).
- **The reference is `c159_<set>`**, all seven sets: the code as committed now, the 2026-09-26 round's
  final candidate.  HPGe did not move in that round: `c159_detx` / `c159_manual` are bit-identical to
  `mergedB_detx` / `mergedB_manual`, the runs of the 2026-09-25 landing.  `mergedB_<set>` is the
  scintillator state before the round.  `ref0924`, `merged0925_*` and `mergedA_*` are history
  (upgrade/Wt4's solver changes had moved the results from the branch's old base: five scintillator
  sets 2936/517/84 x271 -> 2934/506/82 x273).
- **The NoisePlusCurvedPower FWHM form**, sqrt(w0^2 + (w1*(E/661)^s(E))^2) with the log-log slope
  s(E) linear in ln E between its values at the ends of the fitted range, is positive and
  non-decreasing by construction; every class's solves use it.  Against the bounded Bernstein form it
  replaced it measured neutral with heavy churn (2937/508/82 x273 -> 2929/508/82 x274).  It is no
  longer a separate commit: `git show 2fa84c17 | git apply -R` removes it while
  `review/fit-peaks-overhaul` exists.
- **The 2026-09-26 round** added, each a `PeakFitForNuclideConfig` field on for every class but HPGe:
  - **ROIs the solve fits worse than no peaks are confirmed after the observable refit** (item B,
    `final_filter_veto_min_lambda` 16).  An ROI the final filter keeps on its peaks' significance,
    although the solve's model fits it worse than a quadratic continuum-only null, is judged again
    after the refit - but only where that test has power (lambda, the chi2 the peaks' shape leaves
    after projection onto the quadratic, >= 16).  It is kept, rescued on the net spectrum by forward
    selection, or dropped.  Judging the SOLVE's peaks instead dropped real lines the solve mis-sized
    (R500 Yb169_Phantom's 55 keV line, z 315).
  - **Narrow local rescues fit a linear continuum** (`rescue_linear_below_num_fwhm` 4): a quadratic
    absorbs 86 % of one peak's chi2 over 3 FWHM and 97 % over 2.
  - **Zero-activity retry** (item C):
    - x-ray matches can anchor a nuclide, the principal line is ranked by yield, and a nuclide is not
      "zeroed" while the solve fits an ROI of its own lines at >= 2x the principal's z
      (`zero_activity_evidence_anchors`; R500 Pd103_Phantom's 22 keV x-rays);
    - ROIs below 150 keV that the retry dropped are rescued on their own data at the planner's widths
      (`zero_activity_rescue_dropped_rois`; SAM Xe133's x-rays and 81 keV line).
  - **Observable refit**:
    - an ROI whose refit loses every peak goes back to the solve's peaks, which ends delivered peaks
      of negative area (`observable_collapse_restores_solve`; NGH Th228_Unsh 776 keV);
    - wide ROIs are split where the data fall back to the continuum between features
      (`observable_split_valley_fraction` 0.15, `observable_split_valley_min_roi_fwhm` 4.5; R500
      Pu239_Unsh 21-117 keV, Xe133_Sh, Th232_Unsh).
  - **Planner and manual stage**:
    - a Constant continuum is judged from the raw sidebands where netting emptied one
      (`cont_constant_from_raw_sidebands`; R500 U233_Sh 2614 keV);
    - search peaks narrower than half the width model stay out of the manual rel-eff
      (`manual_min_search_fwhm_ratio` 0.5; NGH Ca47_Sh's 1297 keV line);
    - a search peak on one of an ROI's own lines (predicted at >= 10 % of the dominant) no longer
      clips the ROI as an obstacle (`obstacle_own_line_fwhm` 0.5; SAM Cd109_Sh 88 keV).
  - **Strong background lines are delivered as peaks** (`observable_deliver_background_line_z` 5,
    `observable_background_lines_without_background`; the user's decision).
    - Which lines: known background lines (NORM, Cs137/Co60, LaBr3's La138 1436 keV and Ba K x-rays;
      `known_background_line`) seen in the foreground at z >= 5, no wider than 2x the resolution, on
      no model peak, and within a FWHM of a fitted ROI.
    - Without a background spectrum the candidates are foreground search peaks on nine common field
      lines: K40, Th232 239/583/911/2615, Ra226 352/609/1120/1764 keV.
    - A line on a stronger line's backscatter peak, E/(1+2E/511), is skipped (NGH's 662 keV Cs137
      makes a ~188 keV hump).
    - They carry the parent nuclide, the user label "background" and the NORM color.
    - This fixed the most frequent severe review complaint: ROI edges on the flanks of unfitted
      background lines (NGH's Cs137 662 keV, SAM and R500 K40).
  - **Measured and left off**, with the reasons in the fields' docs and `TODO.md`: the item A score
    variants (`refinement_score_*`), the solve width limit and balloon re-solve
    (`rel_eff_auto_fwhm_max_ratio_to_model`, `width_balloon_resolve_ratio`), SAM's channel floor
    (`rel_eff_auto_fwhm_channel_floor_factor`), background lines held undelivered inside ROIs
    (`observable_fixed_background_line_z`), dominant-line admission
    (`admission_dominant_line_required`), the final filter's linear rescues
    (`final_filter_rescue_linear`, `final_filter_rescue_linear_at_search_peaks`),
    `rescue_forward_selection`, `observable_split_valley_by_data`, `observable_inverted_step_as_linear`
    and `observable_skip_at_background_lines`.
  - **Tools**: `invariants.py` checks for overlapping ROIs; `fwhm_balloon.py`; `review/RUBRIC.md`
    has a section for LaBr3 and CZT.
- **The 2026-09-24 round** (74dea9e4) added:
  - **Width model, non-HPGe** (`detail::fit_fwhm_function_robust`):
    - search peaks reaching below the first live channel do not vote (on the GR1 CZT, dead below
      39 keV, a sliced line had set the widths to half the detector's);
    - when no search peak reaches z 6, the class prior sets the scale;
    - when the width model is that prior and the first solve's widths leave [0.6, 1.6] of it, the
      solve is re-run with the widths held (weak GR1 solves had run to 0.13-0.27 of it).
  - **Escape peaks on scintillators** (`escape_peaks_non_hpge`), modelled where the search finds them.
    Tl208's single escape (2103.5 keV; SAM's truth centroid is 2121.5) is now found on SAM (z 22-27)
    and LaBr3 (z 9-15), and Y88's 814 keV double escape on the GR1.
  - **The obstacle fix.**  A group admitted as data-evident on its own search peak no longer treats
    that peak as an obstacle: the nearest search peak within 0.5 FWHM of the line.  It had clipped
    R500 I123_Phantom's Te x-ray ROI above the peak and started Ac225_ShHeavy's 1567 keV ROI on its
    flank.
  - **Three refinement-score variants behind flags, all OFF**, see pending item A:
    `refinement_delivered_segments`, `refinement_whole_range_score`, `refinement_keep_best`.
  - **Tools** moved to `target/peak_fit_improve_ai/`.
- **Numbers** (strong / moderate / weak truth lines found, significant extras; 300 s, Livermore):

  | Set | ref0924 (end of 09-24) | mergedB (landed 09-25) | c159 (now) |
  |---|---|---|---|
  | R500 | 563/103/14 x73 | 563/100/17 x73 | 563/102/16 x86 |
  | NGH | 530/90/17 x35 | 529/91/15 x36 | 537/90/17 x82 |
  | SAM | 683/82/19 x141 | 680/81/20 x136 | 687/81/23 x215 |
  | LaBr3 | 748/142/22 x19 | 745/141/19 x23 | 747/139/20 x63 |
  | CZT | 412/100/12 x3 | 412/95/11 x6 | 412/99/11 x7 |
  | Total | 2936/517/84 x271 | 2929/508/82 x274 | 2946/511/87 x453 |

  - c159's extras include 177 delivered background peaks (R500 14, NGH 44, SAM 78, LaBr3 41), which
    the scorer counts as extras unless the truth lists that line as don't-care; the other extras are
    276 (mergedB 274).
  - HPGe (Detective-X inject, hand fits) was bit-identical through every scintillator round; the
    09-25 landing moved it slightly (Detective-X moderate lines found 555 -> 559, hand fits 470 ->
    469, strong unchanged).  Now Detective-X 991/559/70 x16, hand fits 920/469/14 x22.
  - Visual reviews in the 2026-09-26 round: the baseline R500 review (205 spectra) GOOD 135 / MINOR
    43 / BAD 27.  Changed spectra, better / worse / mixed / same against the previous candidate:
    c137n vs mergedB 38 / 6 / 5 / 20, c155 vs c142 20 / 8 / 8 / 22, background delivery c157 vs c155
    150 / 11 / 14 / 45, c159 vs c157 31 / 2 / 1 / 11.
- **Reference material** in `~/fit_peaks_work/` (`NOTES.txt` there):
  - runs `c159_<set>` (now) and `mergedB_<set>` (as landed on 09-25), for `<set>` in r500, ngh, sam,
    labr, czt, detx, manual; the older `c77_<set>` and `ref0924_<set>` are history;
  - the binaries `bins/eval_c159` (sources in `bins/src_c159`) and `bins/eval_mergedB0`;
  - galleries: `galleries/c159/<detector>/gallery.html` for the five scintillator sets, and
    `galleries/mergedB/{Detective-X,manual}` for HPGe, which has not changed since;
  - review pages: `triage/base/baseline_r500_review.html` (the round's baseline review),
    `triage/c155/c155_vs_c142_review.html`, `triage/c157/c157_vs_c155_review.html`,
    `triage/c159/c159_vs_c157_review.html`, and `triage/c83_wholerange_review/`, the review of the
    whole-range score.

## Where to work next: the plan

**Why this order.**  The 2026-09-26 round's one big win was structural: delivering background lines
as peaks (review 150 better / 11 worse).  Most of its local gates measured neutral and were left off
(the refinement-score variants, the width limits, dominant-line admission, the linear rescues), and
they failed the same way: the fitter has no model of NaI's non-peak structure - the scatter hump, the
backscatter peak, Compton edges, the turn-on knee, the background's own shape.  Low-energy NaI
continua curve on the scale of a peak, so every local test (a quadratic or linear null, a SNIP
baseline, a search-peak gate) calls that structure a peak about as often as it catches a real line.
The refinement lottery (item A), the phantoms and the width balloon (item C) are symptoms of it.
Strong lines are near their ceiling (94-97 % found, SAM 90 %), so expect small moves in the truth
counts; look for fewer phantoms and cleaner fits, which only the visual review measures.  So: first
check that what we have generalizes, then give the fitter structure knowledge, then return to the
lottery.

**Step 0 - ask the user** the open decisions ("Decisions waiting on the user" below).  The scorer
exemption for delivered background peaks matters to every later step: it is a harness-only change
that makes the extras comparable again (177 of c159's 453 are delivered background peaks).

**Step 1 - breadth and the GUI's inputs: measure, do not tune.**
- Run your baseline snapshot on corpora the rules were never tuned on, all in `$FPR_INJ`: the NaI
  D3S and RadEagle (`--det-type=Low`), IdentiFINDER-LaBr3 and LaBr3_1.5x1.5_SNL (`LaBr`), the
  1cm-1cm-1cm CZT cube and CZT_H3D_M400_ORNL_25cm (`CZT`; the H3D resolves far better than the GR1
  the CZT width rules lean on).  Also the R500 at 30 s and 1800 s (`FPR_DWELL`), one set with the
  Denver and Baltimore backgrounds (`FPR_SITE`), and all five sets with `--background=none` (only
  NGH has been measured without a background).  `run_corpus.sh` takes the detector and det-type;
  raise the timeout for 1800 s.
- Build galleries and review ~40 spectra per new set with the rubric (add a detector section to
  `review/RUBRIC.md` where the physics differs).  Tally the defect classes against the five sets'.
  A five-set rule that breaks elsewhere gets fixed before anything new; a new class joins the ranking.
- Cover the GUI's inputs (item H) with harness options, and measure the five sets with each:
  - `--drf NAME`: pass a detector response, as the GUI does when one is loaded; its FWHM curve then
    replaces the class curve as the width prior.  `data/common_drfs.tsv` holds GADRAS responses with
    FWHM for the corpus detectors (IdentiFINDER-R500-NaI, ICx/FLIR IdentiFINDER-NGH, Radseeker-LaBr3,
    Kromek GR1, ...), app-URL encoded.
  - a second fit that is given the first fit's peaks as existing peaks (the user pressing "Fit Source"
    twice): it should reproduce the first fit, not duplicate or shift peaks.
  - `--det-type=auto`: classify as the GUI does (`PeakFitUtils::coarse_det_type`) instead of forcing
    the class.
  The fitter must stay general: a DRF may help, but the no-DRF path remains the one tuned.

**Step 2 - a background-aware observable refit** (get the user's go-ahead on the design first).
- The mechanism: the solve fits the net spectrum (foreground less the live-time-scaled background),
  but `compute_observable_peaks` refits the GROSS foreground with a polynomial or step continuum, so
  the background's structure under an ROI is absorbed by its peaks and continuum.
- Exemplars: R500 U235_Unsh_0400/4000/9000 (the continuum runs ~200 counts below the background's
  50-100 keV hump, and broad "x-ray" peaks take it); LaBr3 Sb124_Unsh 1325/1366 and Tl200_Sh 1408 keV
  (beside La138's 1436/1468 keV feature); NGH W187_Sh (the background's 662 keV line); and the item B
  confirmation, which judges parts holding a background line on the gross spectrum.
- Design: continuum = scaled background + a fitted polynomial.  `PeakContinuum::External` is a fixed
  histogram with no free parameters, so either add an offset type for this, or refit on the net
  spectrum and deliver an External continuum sampled from the scaled background plus the net refit's
  polynomial.  Smooth a sparse background before using it as a shape (a 300 s background has few
  counts per channel).  Settle with the user how it relates to background-line delivery: with the
  background in the continuum, the continuum could draw the background lines instead of delivering
  them as peaks, which is not what a hand fit does.
- Without a background spectrum nothing changes.  Gate it off for HPGe, measure it alone, review.

**Step 3 - physics-predicted structure as phantom evidence.**
- A strong line's own structure sits at energies physics fixes, with no detector response needed:
  the backscatter peak near E/(1+2E/511) and the Compton edge at E - E/(1+2E/511), well apart from the
  line above ~300 keV.  The backscatter rule for background lines (a line on a stronger line's
  backscatter is skipped: NGH's ~188 keV hump from the 662 keV line) was a first instance, and it
  worked.
- Study before coding (Python on `per_peak.tsv` and the trace): of c159's delivered peaks, matched
  and extra, how many sit within ~1 FWHM of a stronger delivered line's backscatter or Compton edge,
  and whether the search peak's width there separates the two.  Already measured: a search peak
  within 0.5 FWHM of a line and at most 1.4x the width model marks a real line ~99 % of the time
  (c135n: R500 442 matched / 6 extra), and humps and shoulders give wide search peaks - but so do
  unresolved x-ray blends, and most extras have no search peak at all.
- The rule to try: at a predicted structure site, a peak the requested sources do not predict
  strongly must show photopeak evidence (a narrow search peak) to be admitted or delivered.
- Exemplars (the c83 review's phantom regions): R500 Np237_Unsh's backscatter, Zr95_Unsh's Compton
  edge, La140_Sh's plateau, Mo99_Sh's scatter hump; SAM's weak background lines on shelves (item E).

**Step 4 - the refinement lottery (item A)**, on top of steps 2 and 3: the whole-range score's
baseline includes the background model, and no ROI earns credit for modelling a predicted structure
site without photopeak evidence.  Re-measure whole-range + keep-best, and review it against the c83
review's BETTER 10 / WORSE 24.

## Pending work: every open item

These are all the open items from `TODO.md`'s sections "fit_peaks_for_nuclides on scintillators",
"fit_peaks_for_nuclides merged onto upgrade/Wt4 - follow-ups" and "fit_peaks_for_nuclides
scintillator round 2026-09-26", which have the full text of each.  They are grouped here, with
exemplars (`det/Spectrum`).  Take them in the order of the plan above; the rest by what the reviews
find most often.

**Decisions waiting on the user** (ask before acting on any of them):
- **Cs137 in the no-background line list?**  Without a background spectrum, delivery is nearly inert
  on NGH (2 lines in 205 spectra, `--background=none`): NGH's strong background line is Cs137 662 keV,
  which is not on the nine-line list.  But an unrequested 662 keV peak could be a real source.
- **Scorer exemption** for delivered peaks labelled "background": they are 177 of c159's 453 extras.
- **SAM-Eagle's low-energy search**: most SAM misses below 100 keV have no search peak at all.
- **A background-aware observable refit** (plan step 2): the solve fits net data, the refit the gross foreground
  with no background model, so background structure under an ROI is absorbed by its peaks (R500
  U235_Unsh_0400: the continuum runs ~200 counts below the background's 50-100 keV hump).  A
  continuum of the live-time-scaled background plus a fitted polynomial would fix it, but delivered
  peaks would then carry a background-shaped continuum.

**A. Refinement acceptance: the "lottery", still the biggest lever.**
- **The mechanism.**  The loop re-plans from each solve and keeps a challenger only if its score
  beats the incumbent's, stopping at the first rejection.  Any change to the plans or the score
  re-decides near-ties, so strong lines move both ways on unrelated spectra.
- **The default common-domain score is biased.**  `PeakFit::chi2_for_region` ignores its channel
  range for continua with an energy range (shadowed locals).  Taking the peaks centred in an ROI pulls
  in the neighbour's copy of a line with the neighbour's continuum.  Both inflate the score of a split
  pair: CZT/Mo99_Sh's 739.5/777.9 keV split scored 6.2 vs 1.5 and lost the z=12.5 778 keV line.
- **Measured, 2026-09-24** (strong/moderate/weak, extras, over the five sets):
  - (a) the corrected common-domain score, `refinement_delivered_segments`: 2915/522/82 x267, noise.
    It favours an ROI fit exactly to the shared range over a wider one holding a line at its edge
    (R500/Br76_Phantom 559 keV), and cannot credit a line a challenger adds (R500/Th228_Unsh 2614 keV).
  - (b) whole-range against the SNIP continuum plus keep-best: 2938/541/86 x282, real lines found
    (NGH/Ca47_Sh 1297, R500/I123_Phantom 159, SAM/Tl200_Sh's seven), but the review of 102 changed
    R500 fits was BETTER 10 / WORSE 24 / MIXED 6 / SAME 62: the SNIP does not follow NaI's scatter
    plateau, backscatter hump, Compton edges and turn-on, so an ROI modelling that structure with
    peaks scores as progress.
  - (c) keep-best with the default score: noise.
- **Measured, 2026-09-26, all left off** (`refinement_score_*`): a narrower SNIP (1.0 and 0.75
  FWHM), its own presmooth, a complexity charge and a power veto bring the c83 WORSE patterns back
  unchanged (R500 La140_Sh's 128 keV plateau phantom, Mo99_Sh's 88-243 keV comb) - the SNIP's
  min-filter rides low on noisy rounded plateau tops at ANY window.  Capping an ROI's credit at its
  gain over a quadratic null also zeroes real lines' credit (I123_Phantom 159 keV z 130,
  Yb169_Phantom 177 keV z 116).  A baseline built from the spectrum alone cannot separate NaI
  structure from peaks.
- **Next ideas**: a structure model (scatter hump, backscatter, Compton edges, turn-on) for the
  whole-range comparison, or health-aware tie-breaking.
- Exemplars to keep in the regression set: CZT/Mo99_Sh, R500/Br76_Phantom, R500/Th228_Unsh,
  SAM/Tl200_Sh, R500/I123_Phantom, R500/I124_Sh, SAM/Uore_Unsh, R500/La140_Sh.

**B. What is left of the worse-than-no-peaks confirmation** (the confirmation itself landed).
- It judges a part holding a delivered background peak on the gross spectrum, where that line can
  carry a weak ROI past the test; a vetoed part's rescue on the net spectrum then drops the
  background peak.  The rescue's continuum is net-level even with a background (pre-existing).
- **No local polynomial null separates NaI structure from peaks.**  A quadratic has no power in an
  ROI of 2-3 FWHM (most scintillator single-line ROIs, every SAM-Eagle ROI); a linear one, or a
  quadratic over the ROI plus a FWHM of sideband, passes real narrow-ROI lines (NGH Pu238_Unsh 153 keV,
  R500 Pd103_Unsh 40 keV) and slope phantoms (R500 Bi213_Phantom 50 keV, NGH Eu154_Unsh 82/94 keV)
  alike, because low-energy NaI continua curve on a peak's scale.  The separating evidence has to
  come from elsewhere: search-peak width (a search peak within 0.5 FWHM and at most 1.4x the model
  marks a real line ~99 % of the time; not yet measured on a full run), the physics model's expected
  amplitude, or a background-aware continuum.
- `compute_roi_chi2_significance` reads `peaks_in_roi[i]->sigma()` with `i` indexing the sorted,
  filtered means, so the dof clustering uses the wrong peak's width.  Shared with HPGe: fixing it
  moves HPGe.

**C. Solve and rel-eff failures.**
- **Solve widths balloon at low energy** (NaI, in both FWHM forms): 46 R500, 67 NGH and 19 LaBr3
  final solves run over 1.6x the planner's model somewhere at or below 300 keV, where the planner's
  model matches the truth widths.  Two fixes were measured and left off: box bounds within 1.4x of
  the planner's curve (R500 13 truth lines gained / 13 lost, SAM lost strong low-energy lines), and a
  re-solve with those bounds when the first solve balloons (32 better / 28 worse).  Narrowing alone
  moves the misfit to flanks, the turn-on and plateau phantoms: the balloon stands in for unresolved
  x-ray blends and the scatter hump.  Next: an upper bound below ~200 keV only, skipped where the
  channel floor sets the widths, or together with a structure model.  Chain example: NGH
  Tl201wTl202_Phantom lost its 167 keV line (widths 27 keV at 26 keV, Tl201 zeroed, calibration slid
  4.5 keV, the rescue measured the lines at the wrong place).
- **RelActCalcAuto's starting curve is not the caller's**: with `starting_fwhm_*` supplied, the seed
  goes through a 3-coefficient response fit to 20 linearly spaced samples, so R500 Co56_Unsh started
  at ~6 keV FWHM at 400 keV where the planner's model says ~40.  Seeding from the caller's curve may
  help every scintillator solve, but it is shared with HPGe: measure it separately.
- **Source x-rays the decay data under-predict** (I123 Te K, In111 Cd K) are simply not data
  evidence.  The right fix is a FloatingPeak at the source's x-ray centroid.  First fix
  RelActCalcAuto's initial estimate, which matches search peaks to lines of any yield (a 100 keV peak
  paired with an I123 line of yield 9e-10).
- **The manual rel-eff ladder forgives lines a form cannot predict**: an excluded peak costs nothing,
  so a 2-point physical fit that excludes a z=53 peak beats the 3-point fits.  NGH/Ca47_Sh's 1297 keV
  line is now found (the narrow-search-peak filter), but the ladder should still charge an excluded
  peak a capped chi2.  Shared with HPGe.
- **First solves moved by upgrade/Wt4's solver changes**: SAM/Am241_Sh lost its 101/123 keV lines
  (z 52/56).  Bisect 20f9d851 / 9e5c18cb / a3d39f4e; these solves are ill conditioned.

**D. Planner and admission.**
- **Moderate groups admitted on a minor line**: R500/W187_Sh's W K x-ray group is predicted at z 53
  with data z 0 at 72 keV and gives a 67 keV phantom.  Requiring a peak at the dominant line was
  measured and left off: a blended x-ray group's dominant line shows no peak of its own either (10
  real lines lost for 6 extras).  The rel-eff misses the shield's x-ray attenuation; another idea is
  needed.
- **The planner joins several groups into one ROI** (R500 Lu177m_Unsh 85-451 keV, Ag110m_Unsh
  524-1014 keV): the comb's peaks fill the valleys, so the valley split cannot cut it.  The planner's
  leak rule has to stop joining the groups.
- **Turn-on above the threshold-ramp floor.**
  - The NaI roll-off knee: R500/Np237_Unsh's 26 keV "peak".
  - Source x-rays on NGH's steep rise: Tc99m_Phantom, U235.
  - The extent-clipped admission gap: a knee-aware admission fixed Np237 but broke Lu177m_Sh's
    neighbour ROI.  Exempting escape groups is the untested next step.
  - SAM's extent estimator returns 0 at 12.5 keV/channel.
  - The NaI/LaBr3 LLD channel: close the admission gap first.
- **CZT/GR1 moderate lines**: the NaI-tuned gates may not suit CZT.
- **Escape peaks the search misses**:
  - SAM/Y88_Unsh's 1343 keV single escape (z=21) sits under a 250 keV-wide search feature.
  - On CZT the missed one is the parent: U233_Unsh and U233_Sh 1592 keV, Th228 1592 keV (z 6.5-9.4),
    where the 2614 keV peak is weak or unfound.  A requested-source line >= 1.6 MeV could stand as
    parent when the search found the escape.  R500/U233_Sh's 2121 keV has the same cause.
- **Phantom low-energy peaks in phantoms**: R500/Yb169_Phantom's 50 keV x-ray is skewed and shifted
  by small-angle scatter.

**E. Observable refit.**
- **The side cap** uses the peak's own FWHM, so a merged x-ray blob (R500/Ir192_Sh's Pb K x-rays as
  one 31 keV peak) trims nothing.  The ruler passed in is the planner's width model, explicitly now
  (it used to be read from a moved-from solution); use it for the edge too.
- **Weak background lines on SAM shelves and Compton edges** (583-599 keV after the 511 keV step,
  911-945 keV on flat stretches) are delivered as phantoms, and no simple test separates them from
  weak real K40 lines: requiring a foreground search peak or a fixed-shape foreground z cost more
  real lines than phantoms (c158).
- **LaBr3 lines next to the La138 1436 keV feature** (Sb124_Unsh 1325/1366, Tl200_Sh 1408 keV): a
  background line ON a model peak is left out (the peak stays gross, InterSpec's convention), so
  Sb124's 1437 keV line slides onto La138's 1468 keV feature and the confirmation, judging the gross
  refit, loses the weak lines beside it.

**F. Detector-specific.**
- **SAM-Eagle's 12.5 keV/channel**: the 58.8/67.7 keV pair shares adjacent channels and a free-width
  refit merges them.  Lowering the solve's channel floor (1.25 channels) barely moved the widths
  (+3 / -4 truth lines), so the floor is not what holds them wide.  Use a fixed-width, fixed-mean
  amplitude fit for sub-channel lines.
- **LaBr3 width scale pulled wide**: La138's internal 1436 keV feature (FWHM ~50 keV, z~30 in every
  spectrum) and the 790 keV beta hump vote k = 1.6-2.0 on 37 of 205 Radseeker spectra (median 1.0).
- **CZT single-ROI energy-cal slide**: CZT/Pd103_Unsh's one ROI on the threshold-sliced 39.8 keV
  line let the energy offset slide 8.7 keV, below the live channels.  With one ROI the offset should
  stay fixed.
- **GR1 physics, for the record**: the GR1 is dead below 39.3 keV.  Pd103's 20-23 keV Rh x-rays have
  no data; its 39.8 keV line is cut in half; no fit is the right answer there.
- **Scorer**: truth lines centred below the first live channel should be dontcare.  GR1/I125_Unsh's
  28.6 keV (z=244) is scored as missed.
- **Truth lines InterSpec cannot own**: NGH/Np237_Sh's "405.5 keV" line (z 130, 58 % of the 312 keV
  line; the decay data give ~1 %), NGH/Xe133_Sh's 302.9 / 383.9 keV lines (Ba133's energies).  Worth
  checking the inject truth.

**G. Code hygiene.**
- `chi2_for_region`'s shadowing: the HPGe path still uses it.
- `RelActCalcAuto::Options` XML fields added without bumping `sm_xmlSerializationVersion` (now also
  `fwhm_max_ratio_to_start` and `fwhm_channel_floor_factor`).
- `recovered_source_anchors` is a const empty vector, so every "recovered anchor" check in the
  refinement is a no-op: wire it or delete the checks.
- `PeakFit::fit_amp_and_offset_imp` returns the chi2 with the continuum clamped at zero.  Forward
  selection is designed on that value (the unclamped chi2 rewards negative continua), but the
  backward elimination and other callers should be audited.
- **Speed**: after the merge R500 Uore_Unsh took 157 -> 284 s and As72_Unsh 5 -> 23 s; profile them.

**H. Breadth, not yet measured.**
- Only five detectors at 300 s / Livermore have been used.  The inject corpus has ~30 more
  scintillators (`ls $FPR_INJ`), plus 30 s and 1800 s dwells and Baltimore/Denver backgrounds
  (`FPR_DWELL`, `FPR_SITE`).  The no-background path was measured on NGH only.
- The CZT width rules lean on the class curve, which IS the GR1's (`czt_fwhm_fcn`), so the
  high-resolution CZT_H3D_M400_ORNL_25cm is the test of their generality.
- Try at least one more NaI (D3S, RadEagle), one more LaBr3 (IdentiFINDER-LaBr3) and a CZT cube.
- **The GUI's inputs**: the Reference Photopeaks tab's "Fit Source" runs this fitter with the same
  configuration, but with inputs the corpus never used: the loaded DRF (whose FWHM curve replaces the
  class curve as the width prior), the Peak Manager's existing peaks (replaced where a fitted ROI
  covers them), the automated search run with that DRF, and the detector class InterSpec guesses
  from the file.  A corpus variant with a DRF and with pre-existing peaks would cover them.

**The round's baseline R500 review** (`triage/base/baseline_r500_review.html`): severity >= 2
classes CONT 29, PHANTOM 12, MISFIT 11, EMPTY 9, MISSED 9, COMB 9, MERGE 8.  CONT's most frequent
cause, ROI edges on the flanks of unfitted background lines, has been fixed since; review c159 again
before choosing.

## Measurement kit

Everything is in `target/peak_fit_improve_ai/` (`README.md` there is the index): the harness in
`harness/`, the scripts in `tools/` (they source `tools/env.sh`), the rubrics in `review/`.

- **Setup, once.**
  - `python3 -m venv ~/fit_peaks_work/venv && ~/fit_peaks_work/venv/bin/pip install numpy matplotlib`,
    then `export FPR_PY=~/fit_peaks_work/venv/bin/python`.  The system python3 has no matplotlib.
  - The work area defaults to `~/fit_peaks_work`.
- **Build:** `cd target/testing/build_ninja && ninja fit_peaks_corpus_eval test_fitPeaksForSources test_RelActCalcAuto_ProfileApi test_fitPeaksCorpusScore`.
- **Snapshot** every build you will measure: `tools/snapshot.sh cNN`, continuing from c160.  Measure
  only `~/fit_peaks_work/bins/eval_cNN`: other Claude sessions share `build_ninja` and rebuild it
  mid-run.
- **Corpora:**
  - `$FPR_INJ/<detector>/<site>/<dwell>_seconds` holds the GADRAS-injected spectra with truth.  The
    five sets: IdentiFINDER-R500-NaI, IdentiFINDER-NGH, SAM-Eagle-NaI-3x3 (`--det-type=Low`),
    Radseeker-LaBr3 (`LaBr`), Kromek-GR1-CZT (`CZT`).
  - **The harness defaults to High**, so a scintillator run without `--det-type` silently evaluates the
    HPGe configuration.  The scripts pass it, plus `--background=file --no-structure
    --weight=min_scored_energy=20` and a 1200 s timeout (U and Pu spectra take 2-5 minutes).
  - `run_hand.sh` runs the 17 R500 hand fits (geometry against a human).
- **Runs:**
  - `quick.sh BIN OUT DET TYPE a,b,c` renders images for a subset; set `FPR_CMP=<old run>` for
    red-new over blue-old.
  - `full.sh BIN TAG [--set f=v]` runs the five sets plus the HPGe sets (1.5-2 h).  Add
    `FPR_GUARD_REF=c159` (or your baseline's tag) for the HPGe bit-identity check.
  - `--background=none` (the last `--background` argument wins) runs the no-background path.
  - `gallery.sh` builds HTML galleries.
  - `--problems a,b,c`, `--set field=value` and `--debug ID` (one problem, serially, with the fitter's
    own trace).
- **Outputs:**
  - `per_peak.tsv`: verdicts for every truth and fitted peak.
  - `per_problem.tsv`.
  - `roi_plan_trace.txt`: every line group admitted or rejected and why, share decisions, continuum
    choices, and one line per solve and refinement challenger with its score.
  - `plot_data/*.json`.
- **Tools:**
  - `cmp_runs.py OLD NEW`: gained/lost truth lines; the main A/B.
  - `run_table.py PREFIX...`: the per-set table above.
  - `extras_cmp.py OLD NEW`: read it whenever lines are gained.
  - `miss_why.py`, `gate_rej.py`, `refine_margins.py`, `worse_delivered.py`, `invariants.py`,
    `bitident.py`, `fwhm_balloon.py` (solve widths against the planner's model).
  - `changed.py OLD NEW`: which spectra changed; the ids are the lines starting with a letter.
  - `render_all.py` / `overview.py` (`--cmp`, `--ref`) and `plot_fit_peaks_rois.py` for images.

## Visual review with agents

- **Render** the changed spectra (`changed.py`, then `render_all.py RUN OUT --cmp OLD` or
  `overview.py`).  Look at the exemplars yourself with the Read tool.
- **For volume**, split the ids into batches of ~20 and launch one general-purpose agent per batch,
  5-8 in parallel.  A batch costs 250-300k tokens.
  - Rubric: `review/RUBRIC.md`, written for NaI; it has the iodine-escape physics notes without which
    reviewers mis-call real peaks, and a last section for LaBr3 and CZT.
  - Plus `review/RUBRIC_BOTH.md` for new vs old.
- **Tally** with `triage_sum.py FILE...` and read the WORSE notes.  `review_page.py` makes an HTML
  page; `rr_compare.py` compares a re-review with the review before it.
- **Check "phantom" or "junk" claims against `per_peak.tsv`** before acting.  Reviewers call real
  faint lines phantoms.
- Agent prompt (substitute the paths, `<N>`, `<NEW>`, `<OLD>`):

  ```
  You are reviewing automated gamma-spectrum peak fits (<detector type> detector) by LOOKING at images,
  as an expert spectroscopist would.  Instructions: <repo>/target/peak_fit_improve_ai/review.
  1. Read `RUBRIC.md` there first and view its three calibration images.
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
   - Snapshot the current tree and run `full.sh` with `FPR_GUARD_REF=c159
     FPR_GUARD_SETS="r500 ngh sam labr czt detx manual"`.  If nothing on upgrade/Wt4 has touched the
     fitter since 2026-09-27, every set is bit-identical to `c159`; if not, find out what changed,
     and use your baseline as the reference from then on.
   - The five sets' galleries are in `galleries/c159/`, and R500 had a full review at the start of
     the last round; review a sample of each set at c159, and the new sets of the plan's step 1.
2. **Triage** the review notes and `miss_why.py` reasons into classes; follow the plan above, and
   let what the reviews show re-rank it.
3. **Trace** 2-3 exemplars through `roi_plan_trace.txt`, then `--debug`, until you can name the code
   path.  Temporary hooks are fine; remove them before full runs (grep the diff for `getenv`,
   `TEMPORARY`).
4. **Fix the mechanism, gated.**
   - A `PeakFitForNuclideConfig` field, OFF in the header, ON in `s_default_non_hpge_config` (which
     covers NaI, LaBr3 and CZT).
   - Registered with `FPN_*_FIELD`.
   - The *why* and the exemplar in the comment.
   - Shared code (`PeakFitLM`, `RelActCalcAuto`, `PeakFit`) changes are opt-in options, off by default.
5. **Quick check** the exemplars plus the regression set with `quick.sh`, and look at the images.
6. **Full check**:
   - `full.sh` on every set and `invariants.py`;
   - HPGe bit-identity;
   - both unit suites (add a unit test for a new mechanism);
   - re-review the changed spectra, and compare with `rr_compare.py`.
7. **Record** in the memory note `project_fit_peaks_nai_state.md` and in `TODO.md`, then take the
   next class.

## Rules learned the hard way

- **Noise.**  ±2-4 strong lines per set between candidates is noise (the refinement lottery).  Read
  the gained/lost lists; when a strong line flips, re-run it on both binaries with `--debug` and find
  the step that changed.
- **Determinism**: a difference between runs of the same binary is a code difference or a timeout.
- **Look.**  More truth lines is not the goal.  The whole-range score (+51 lines, twice as often
  worse) is the example.
- **A refactor must leave the default path bit-for-bit alone.**  Re-run 3-5 problems and compare the
  exact trace numbers (refinement scores) with the reference run before a long run.
  `solution_chi2_over_segments`' `delivered_model` argument already chose the significance test.
  Reusing it for a second purpose silently changed every NaI refinement for two full runs.
- **Measure a change alone before stacking it**, and run every set after a planner gate change.  A
  data test tuned on R500 (2.9 keV/channel) broke SAM (12.5 keV/channel) and NGH.
- **Planner geometry changes re-decide the refinement passes.**  Fix presentation (ROI widths, joins)
  downstream, in the observable refit.
- **A rule that deletes a phantom ROI can break its neighbour**, whose boundary it had held
  (Lu177m_Sh).
- **A fix can remove an accidental benefit.**  The obstacle bug had been cutting I123_Phantom's
  x-rays out of a solve that could not fit them together with 159 keV.
- **A window rule must not reach the next feature.**  The first obstacle fix removed any search peak
  within a FWHM, took Tl201woTl202_Phantom's 61 keV neighbour with it, and lost a 40.7 keV line.
- **Truth quirks** (GADRAS):
  - Truth includes iodine and pair-production escapes.
  - Energies can be blend centroids.
  - Some sources lack lines the decay physics predicts (Mo99 without Tc99m).
  - Some truth lines disagree with the decay data.
  - InterSpec labels lines by the ultimate parent; never "correct" that.
- **CZT's discriminator**: check where the first live channel is before trusting a low-energy peak.
- **What is not delivered cannot be drawn.**  The solve fits net data, the observable refit the gross
  foreground.  Background lines held in the model but not delivered left the drawn fit below the
  data, and the review called that worse (6 better / 39 worse) although the areas improved.
- **Local nulls.**  A quadratic continuum absorbs 86-97 % of one peak at 2-3 FWHM, so a local
  peaks-vs-null test there has no power; a linear one has power but calls NaI's curved low-energy
  continua (knee, scatter hump, Compton shoulder) peaks.
- **`PeakFit::fit_amp_and_offset_imp`'s chi2 holds the continuum at or above zero.**  Backward
  elimination from every line then keeps nothing when one line needs a negative continuum; scoring
  the unclamped chi2 instead rewards negative continua.  Forward selection on the clamped value works.
- **Review a stacked candidate step by step.**  c155's review split 4/0/1, 12/6/6 and 5/2/2 (better /
  worse / mixed) by step showed which change caused the losses.
- **The observable refit runs ROIs on a thread pool**: `record_roi_plan_trace` from inside it is lost
  (thread_local); record into the notes collected after the pool joins.

## Code map

- `src/FitPeaksForNuclides.cpp`:
  - `plan_rois_impl` (the planner): groups, admission, data-evident groups and obstacles, sharing,
    span cap, extents, turn-on floor, continuum type.
  - `estimate_initial_rois_using_relactmanual`: the manual-stage first plan.
  - `fit_peaks_for_nuclide_relactauto`: the first solve and the width-runaway re-solve; retries
    (no-ecal, physical model, zero-activity); the refinement loop (`score_pair`, anchor guards,
    keep-best); the final filter and rescue; observable peaks.
  - `compute_observable_peaks`; before it, `add_background_lines_to_model`, `known_background_line`
    and the common-line tables (`sk_common_background_lines`).
  - `rescue_roi_locally` and `quadratic_null_power`: the local rescues and the final filter's
    confirmation.
  - `detail::fit_fwhm_function_robust`: the width model.
  - `solution_chi2_over_segments` / `delivered_spectrum_chi2`: the refinement scores.
  - `add_escape_peak_floating_peaks_if_appropriate`.
- `InterSpec/FitPeaksForNuclides.h`: `PeakFitForNuclideConfig`, every field documented with the
  measurements that set it; `default_config()` and `s_default_non_hpge_config` are in the .cpp.
- `src/RelActCalcAuto.cpp`: the joint solve.  `src/PeakFitLM.cpp`: the refit.
  `src/PeakFit.cpp`: search and `chi2_for_region`.

## Where to read first

- `target/peak_fit_improve_ai/README.md`, and `TODO.md`'s scintillator section.
- Memory notes:
  - `project_fit_peaks_nai_state.md`: every round, mechanism and exemplar.
  - `feedback_fit_peaks_defect_driven.md`, `feedback_fit_peaks_generality.md`.
  - `project_fit_peaks_corpus_harness.md`, `project_fit_peaks_escape_peaks.md`.

## Reporting

After each candidate, report to me:
- a numbers table per set against the previous candidate and the baseline;
- the visual review counts;
- each change, the mechanism it fixed, and its exemplar;
- what got worse, and why;
- the open items.

When a gallery or review page is ready, `open` it for me.  Keep going until the defect classes run
out or I stop you.
