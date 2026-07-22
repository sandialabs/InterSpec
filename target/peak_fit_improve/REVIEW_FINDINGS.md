# fit_peaks_for_nuclides / NuclideConfig GA — Review Findings

Source: multi-agent read-only review (2026-06-03) of `src/FitPeaksForNuclides.cpp`,
`src/RelActCalcAuto.cpp`, `src/RelActCalcManual.cpp`, and `target/peak_fit_improve/NuclideConfig_GA.*`
+ `PeakFitImprove.cpp`.  58 candidate issues → 44 confirmed after adversarial verification.

This is a working punch-list.  Line numbers are approximate (from the review snapshot) and may have
drifted with edits — search by symbol/description.  Status is hand-maintained.

---

# 2026-07-19 — Build made STATIC + eval-nondeterminism hunt (OPEN)

## Build: libInterSpec is now STATIC for the dev tools  [LANDED]
`PeakFitImprove` was a thin loader linking `@rpath/InterSpec.dylib`; the dylib was overwritten on every
build, so re-running any `*_verify` binary silently ran the CURRENT code, and hashing the executable
was meaningless.  Fix: `target/peak_fit_improve/CMakeLists.txt` now sets `USE_BATCH_CLI_TOOLS OFF` and
`USE_MACOS_QUICK_LOOK OFF` → `INTERSPEC_LIB_TYPE STATIC` (top-level CMakeLists:164).  Now the exe embeds
the code (50.9 MB), has no dylib dep, its hash tracks source changes, and builds are reproducible
(verified: flip a constexpr → hash changes; revert → hash returns).  Exposed + fixed an ODR bug:
`LowHighResClassifySettings::print/to_json` were defined twice (src/PeakFitSpecImp.cpp AND, via the
`HighResClassifySettings` alias, ClassifyDetType_GA.cpp) — dup removed from the dev-tool file.
**A/B rule:** rebuild + compare `LibInterSpec/InterSpec.a` (or the exe) per arm; run arms sequentially.

## Eval nondeterminism  [RESOLVED 2026-07-19]
The NaI eval was NOT bit-reproducible run-to-run even on the FROZEN STATIC binary (6 runs over
subsample-4 gave cost 375–478, 3–5 failures; intermittent on `Uore_Unsh|*`, `As72_Phantom|Mirion`).
Two sources; both fixed.
- **Source 1 (was already fixed):** R6 interferer nuclide list iterated from a `std::set<const Nuclide*>`
  in pointer-address order → sorted by `(Z,A,isomer)` (`FitPeaksForNuclides.cpp` ~8185).
- **Source 2 (the real one) — WALL-CLOCK SOLVER TIME CAPS.** The "ASLR/uninit" prior was wrong for the
  eval: it calls `fit_peaks_for_nuclides` with EMPTY `user_peaks`, so every pointer-keyed reorder path
  (incl. `roi_groups` at `FitPeaksForNuclides.cpp:7471`, gated by `!user_peaks.empty()`) is skipped — the
  setup is deterministic.  The cause was that every ceres solve set a *wall-clock*
  `max_solver_time_in_seconds`; under variable CPU load the solve stops after a LOAD-DEPENDENT iteration
  count → different result each run.  Env-gated `[SOLVETIME]` instrumentation (since removed) proved it:
  even SINGLE-THREADED the RelActCalcAuto **pre-fit 20 s cap** tripped ~10×/run at 79–252 iters for the
  same 20 s; localized to `Uore_Unsh|GR135Plus` (serial cost 0.869 vs 4.278 run-to-run).
  **FIX:** raised all wall-clock caps to 1200 s (pathological backstop only) so the deterministic
  `max_num_iterations` governs — `RelActCalcAuto.cpp` main(120→1200)+pre-fit(20→1200), `PeakFitLM.cpp` ×2
  (→1200), `RelActCalcManual.cpp` (60→1200).  **Verified bit-identical:** frozen static binary,
  subsample-40, serial ×2 AND multi-thread ×2 → total_cost 111.107, identical per-spectrum TSVs, 0 caps
  hit (slowest solve now an 81 s pre-fit at its 300-*iteration* cap).  Changes uncommitted.
- Note: `--number-threads 1` DOES serialize the eval; each fit is internally single-threaded
  (`set_max_solve_threads(1)`), so ceres FP-partition order was never a factor.

---

# 2026-07-19 — Change 3 (found-peak seeding): KEEP (seeding HELPS). Earlier "DROP" verdict RETRACTED.

All changes in the working tree, NOT committed.

## ⚠️ RETRACTION of the DROP verdict below
The evidence in the "## Decision / ## Evidence" subsections that follow is **WRONG** and is retracted.
Root cause: **`PeakFitImprove` is a thin loader; the fit code lives in `InterSpec.dylib`, which is
overwritten on every rebuild.** Hashing the *executable* is meaningless (it never changes — `seeded_verify`,
`noseed_verify`, and every fresh build share one executable SHA); the DYLIB is the artifact that carries
the code. So the original "seeded vs noseed" A/B never controlled which code actually ran — both arms
loaded whatever dylib was current, which was **seeding-ON in both** (that is why 356 and 378 were close;
the +22.5 was pure run-to-run jitter). The distribution table below is two runs of the *same* code.

## Corrected finding — seeding is a large NET WIN
Controlled experiment (same executable; dylib rebuilt seeding-OFF vs seeding-ON **sequentially** so each
arm truly loads its own code; subsample-4 / 346 spectra; MT ×3 each; per-spectrum TSV):

| seeding | mechanical failures (×3 runs) | total_cost (×3) |
|---|---|---|
| **OFF** | 6, 7, 6 | 539.5, 575.2, 537.2  (~550) |
| **ON**  | 1, 3, 2 | 366.7, 376.5, 373.2  (~372) |

Seeding ON removes ~5 **"L-M solving failed — no convergence"** failures (As72_Unsh, Tl201wTl202, Co56,
As72_Phantom…) — the found+matched-peak spectra that need the seeded tight ROIs to converge. So the
2026-07-19 morning act of *gating seeding off* was itself the 356→~550 regression. **FIX applied: the
`sm_enable_found_peak_seeding` gate was removed entirely; seeding is unconditional again** (dylib back to
the seeding-ON variant `1d330ee5…`). Change 1 (output-gate OR `max_peak_significance`) and Change 2
(valid-empty) stay.

**Build-hygiene rule (critical):** to A/B any code change, freeze/compare **`InterSpec.dylib`** per arm —
rebuild the dylib for each arm and run them *sequentially* (a parallel run loads whichever dylib is on
disk at `dlopen` time). Never rely on separate-named executables; they all load the one live dylib.

## Open: a residual ~1% run-to-run nondeterminism (real bug, under investigation)
Even seeding-ON varies 1/3/2 failures run-to-run (and ~1% cost) — the L-M problem setup is not
deterministic. Symptom (intermittent no-convergence) points at order-dependent problem construction
(ASLR-ordered pointer-container iteration) or an uninitialised value — NOT threading/FP per the code
owner. Characterization (single- vs multi-thread, failing-set diff) in progress.

---

<details><summary>RETRACTED original DROP evidence (kept for the record — do not trust)</summary>

## Evidence (NaI eval, 346-spectrum holdout: subsample 4, frac 0.2, role eval, seed 20260703)
Seeded vs noseed, compared as **order-independent distributions** over the identical spectrum set
(see the methodology note — a per-spectrum join is impossible with the current TSV):

| metric | noseed | seeded | Δ |
|---|---|---|---|
| aggregate fg cost | **356.3** | **378.8** | **+22.5 (~+6 %, worse)** |
| mean miss_fraction | 0.110 | 0.111 | ~0 |
| total-miss spectra (miss=1.0) | 30 | 30 | **0** |
| zero-observable spectra | 55 | 55 | **0** |
| non-Success | 1 | 1 | **0** |
| Σ n_fit (all peaks fit) | 2523 | 2501 | **−22 (fewer, not more)** |

The n_fit histogram is bin-for-bin nearly identical.  So seeding does **not** rescue empties (30→30),
does **not** proliferate spurious peaks (Σ n_fit falls), and does **not** cause net collapses
(zero-obs 55→55).  Its only population effect is the **+22.5 cost** — diffuse perturbation of the
fitted AREAS of already-fit peaks, via the tight-ROI → rel-eff/activity coupling.  For a change whose
whole purpose was "recover missed peaks / make them all fit," a population result of *empties
unchanged, peaks-fit slightly down, cost up 6 %* fails its own goal.

## Methodology note (a real error I made and corrected)
My first pass **joined seeded↔noseed on `(source,detector)`** and reported a dramatic "60 worse / 59
better, 10 collapses, 17 regressions vs 16 recoveries" churn.  **That was wrong** — there are only 285
unique `(source,detector)` keys across 346 rows (up to 3 spectra share a key: different injects/
shielding), and the two parallel runs write rows in **different completion order**, so the join was a
**cross-product**, not a 1:1 pairing.  The "recoveries" and "collapses" (e.g. `Co56|Mirion 3→0`,
`Ta182|IdentiFINDER 0→2`) were artifacts of pairing a good inject's row with a bad inject's row of the
same key.  The TSV carries **no unique spectrum id**, so per-spectrum pairing is impossible from this
data; only order-independent aggregates/distributions are valid.  **Follow-up:** add the source file
path (or inject index) as a TSV column so future A/Bs can pair per-spectrum.

## Approved follow-up — GATED rescue-seed (design chosen 2026-07-19, NOT yet implemented)
Keep the genuine from-empty recoveries (Sc46/Zr89 — found+matched peaks the keep-gate dropped, leaving
NO ROIs) without perturbing the ~290 already-working spectra that produced the entire +22.5.  Chosen
gate (user-approved): **seed only when the un-seeded ROI set is EMPTY**, checked in-pipeline at the two
existing seed sites — change the guard from `if( sm_enable_found_peak_seeding )` to:
  - manual ROI-estimation site (`estimate_initial_rois_using_relactmanual`): `if( initial_rois.empty() )`
  - auto re-cluster site (before the existing `if( refined_rois.empty() )`): `if( refined_rois.empty() )`
Why safe: working fits (non-empty ROIs) skip seeding entirely → the +22.5 cannot recur; empty fits can
only gain peaks or stay empty (nothing to perturb).  Change 2 (valid-empty) is preserved: a genuinely
empty spectrum has no `peaks_matched`, so the seed helper adds nothing and it correctly stays empty.
**Acceptance bar:** a clean same-batch A/B must show total-miss count < 30 (real recovery) AND
aggregate ≤ noseed (no perturbation) — strictly better than both seeding-on and seeding-off.
Rejected alternatives: retry-on-empty-result (double-fit, more code) and per-source-empty (only matters
for multi-source app fits; ≈ empty-ROI in the single-source eval).

## Reproducibility caveat (2026-07-19) — the +22.5 is softer than first stated
The fast NuclideConfigEval scorer is **not bit-reproducible** for NaI: two identical runs of one frozen
binary gave 164.9 vs 163.4 (89 spectra, ~1% jitter).  On 346 spectra that is ≈±6-7 on a seeded-vs-
noseed pair, so the raw +22.5 is ~3× the noise band, not decisive.  Worse, the **original 356.3
(1.03/spectrum) baseline is unreproducible** — the current build and the frozen `*_verify` binaries all
cluster at ~1.8/spectrum on the identical set; ruled out run-to-run jitter (~1%) and the populated
user-DB (~6%, moves 1.85→1.74/spec).  It was a transient binary/config state from the ~01:00 build-and-
revert window.  **Consequences:** (1) the Change-3 verdict rests on the *structural* evidence (total-
miss 30=30, zero-obs 55=55, Σ n_fit slightly down — seeding does no structural good) + the fact that the
current seeding-off build (159.4) is not worse than noseed_verify (164.9), NOT on the absolute +22.5;
"seeding shows **no reliable benefit**" → gate it off.  (2) Every future accuracy A/B (gated-rescue,
R6 inc 2-3) must rebuild ALL arms fresh, run them back-to-back in identical CWDs, compare only within
that batch, and prefer `--number-threads 1` to kill the jitter.  See memory
[[nuclideconfigeval-reproducibility]].

**Clean same-batch A/B CONFIRMED it (2026-07-19, 346 spectra, identical fresh CWDs, 10 threads):**
seeded_verify **581.9** vs noseed_verify **561.0** = **+20.9** (reproduces the original +22.5 almost
exactly - the seeding penalty is a robust RELATIVE effect even though absolutes drifted 356->561);
newgate (current build, seeding gated off) **555.4** is noseed within jitter -> gate is correct, not a
regression.  **Current trustworthy NaI baseline (subsample 4 / 346 / `--config-genes default` / current
build): ~555.**  The absolute 356->561 drift of the *same frozen binary* is NOT the data (inject +
static data unchanged since 07-17) and only ~6% DB - root cause unresolved (most likely the original
356 log came from a different transient binary than the `*_verify` snapshots); flagged, non-blocking,
since every decision rests on the reproducible relative delta, not the absolute.

</details>

---

# 2026-07-18 (cont.) — FailedToSetupProblem / low-count-NaI "silent empty": the significance funnel

All changes in the working tree, NOT committed.  This continues the low-res diagnosis below.

## The disease: three significance judges that disagree
A source gamma passes through THREE differently-formulated significance tests; the two that gate the
pipeline both systematically dilute per-peak detectability and discard evidence the auto-search
already found:
1. Auto peak-search + source match (`:~6495`) — data peak-finder on FOUND data peaks.
2. **Keep-gate**, `cluster_gammas_to_rois` (`:~4737`): `z = S/sqrt(S+B)`, **S = PREDICTED**
   (`numPerSecond·rel_eff`, so it inherits activity-estimate error), B over the WIDE `±roi_core_num_fwhm`
   core.  Config threshold 5.27/6.39 (a known-bad, objective-inverted GA run; struct defaults 2.0/3.0).
3. **Output-gate**, `compute_roi_chi2_significance` (`:~2961`): Wilks `equivalent_z`, **dof = #peaks in
   ROI** (`:3123` is the sole gate).  A real 6-7σ line inside a 95-peak NORM-chain ROI (Th232) is
   referred to a chi^2_95 null and scores NEGATIVE.
The HONEST per-peak ±1-FWHM Poisson significance (matches the truth's `NSigmaOverBkg`) is already
computed as `max_peak_significance` (`:3113-3120`) and was thrown away (used only for debug).

Measured funnel (RadEye/Livermore/300 s, det=Low), all via the fixed `NuclideFitDebug` harness:
| case | truth σ | auto found+matched | keep-gate z (thr) | output equiv_z | max_peak_sig | observable (before) |
|---|---|---|---|---|---|---|
| Sc46 | 6.2/4.1 | yes @889 | 4.2/4.8 (<5.27) | 1.40 | 2.60 | 0 |
| Co58 | 11 @810 | yes @810 | manual 5.8✓, auto 4.8 (<6.39) | 4.27 | **5.07** | 0 |
| Zr89 | 10 @909 | yes @909 | 4.5 (<5.27) | 3.91 | 4.88 | 0 |
| Th232| 6-7 | yes @86/238/587 | **PASS** 6.9/5.9 | **-1.8/-1.9** | 6.58/7.22 | 0 |

## What LANDED / is happening
- **Scorer split (Step 1 #1, landed):** the objective (`PeakFitImprove.cpp:1006` scorer) scores fit
  QUALITY only; a failed/empty fit is a total-miss via a shared `accuracy_cost` helper (no magic
  penalty).  Mechanical failures are a SEPARATE reliability count, split into real failures
  (miss_frac>0 — had def-wanted peaks) vs legitimate empties (nothing to fit).  Per-spectrum TSV
  (`PEAKFIT_PERSPEC_TSV`) emits status enum, n_observable, n_fit, error.  Full-holdout: total_cost
  7990→3273; 35 real mechanical failures + 16 legit empties.
- **Change 1 (LANDED, confirmed):** output gate at `:3123` = `(equivalent_z >= thr) ||
  (max_peak_significance >= thr)`.  Recovers Co58 (5.07) and the Th232 output-gate class — the class
  the keep-gate sweep provably CANNOT fix (Th232 passes the keep-gate).
- **Keep-gate threshold experiment — ABANDONED.**  A narrow low-count sweep said "lower z to 2.0";
  the broad A/B (346 spectra, all sources) REVERSED it: uniform z=2.0 is +18.7% worse, spawning
  spurious peaks on busy chains (Th232 0→13 fit peaks) the output gate does not fully filter.  The
  `(f·B)²` continuum-systematic term over-suppressed real low-count peaks.  Keep-gate is back to
  clean `z = S/sqrt(S+B)`; scaffolding reverted.
- **Change 2 (LANDED):** make a genuinely-empty ROI set a valid `Success`+0-peaks (not
  `FailedToSetupProblem`); DELETE the keep-gate rescue (a hollow Success-with-0-peaks that also
  triggered an activity-collapse loop: wide ROI → activity 501→13 → auto re-cluster sees all sub-floor
  → ROIs vanish → areas ~half truth → max_peak_sig too low).
- **Change 3 (found-peak seeding) — TRIED, then DROPPED (gated off 2026-07-19).**  Seeding a tight
  (±1.5 FWHM) ROI for every found+matched auto-search peak looked SELF-LIMITING in spot checks
  (Sc46/Zr89/Co58 recovered), but the full A/B verdict is a net **+22.5** aggregate regression with
  **no** population-level benefit — see the 2026-07-19 section at the top.  Now gated behind
  `sm_enable_found_peak_seeding = false`; helper kept for a future empty-baseline-only rescue.

**Verification harness:** `NuclideFitDebug --file <det/city/300_seconds/<src>_Sh.pcf> --det Low --src
<nuc> --debug` (fixed to fit compact meass[0] + subtract meass[2]).  Eval: `NuclideConfigEval
--det-type NaI` prints a mechanical-failure count; drop `--eval-html` (a ~1hr chart render).

---

# 2026-07-18 — NaI crash fix, gross-scoring review + R1.1, low-res diagnosis

All changes below are in the working tree, NOT committed.

## LANDED — `findROILimit` heap out-of-bounds (SIGBUS) on low-resolution spectra  [was blocking all NaI eval]
`NuclideConfigEval --det-type NaI` crashed partway through scoring (SIGBUS / EXC_BAD_ACCESS).  It was
**mis-diagnosed first as a stack overflow** (a `BigStackThreadPool` was added, then reverted).  The
crash reports show the fault address is **~31 GB from the crashing thread's stack pointer** — a wild
**heap read**, not a stack guard-page hit.  Root cause: the **second scan loop in
`findROILimit` (`src/PeakDef.cpp`)** terminates only on `channel != lastchannel && channel !=
meanchannel`, with **no `[0, nchannel)` bound**; on a degenerate low-res spectrum where the scan
direction points away from both targets, `channel` (a signed `int`) runs to `INT_MAX` and
`contents[channel]` reads far out of bounds.  (Tell-tales: the very next line bounds-checks
`channel+direction` but not `channel`; and there is an `assert(channel < num_gamma_channels()+100)`
compiled out in Release.)  **Fix:** hard-bound both scan loops with `channel >= 0 && channel <
nchannel` — changes behaviour only for the currently-UB runaway case.  A first patch that clamped
`startchannel` was reverted as a no-op (`Measurement::find_gamma_channel` already clamps its result).
**Verified:** a full 1464-spectrum NaI scoring pass now completes with no crash; `test_fitPeaksForSources`
is green (no regressions).  This is a genuine engine bug that also affects the app (fits on worker
threads), not just the eval.

## LANDED — Gross-scoring metric: review + R1.1 (truth-area relative-error floor)
Full metric review in `GrossScoring_Review.md` (Opus, read-only).  Headline (confirms the user's
concern): truth `ExpectedPhotopeakInfo::peak_area` is the **analytic, noise-free, full-energy area
GADRAS predicts** — not a fit to the data — and can exceed the counts actually in the ROI (skew; NaI
/ low-energy x-ray efficiency & resolution differences; `create_expected_photopeak` already *warns*
`peak_area > total_area`).  The area term `|A_fit − A_truth| / sqrt(A_truth)` then **over-penalizes a
correct fit** (a 75k-count peak hits the cap of 20 at only ~7.3% truth bias) and lets bright biased
peaks dominate the summed objective.
- **R1.1 (landed):** relative-error floor `sigma = sqrt(max(1,A) + (f_rel·A)^2)`, `f_rel = 0.05`, in
  BOTH `FinalFit_GA::calculate_final_fit_score` and `eval_final_peak_fit`.  The floor only bites when
  `f_rel·A ≳ sqrt(A)` (large-area/biased peaks); dim clean peaks unchanged.  `f_rel` flat for now;
  energy-/det-type-aware is the planned reassess.
- **Deferred (from the review):** R1.2 cap truth at realized `data_area_in_roi`; R2 failure-cliff
  rescale + per-spectrum normalization + two-axis reporting; R3 smooth the 30/50 keV and `tau_def`
  thresholds into continuous reliability weights.  A clean R1.1 A/B (with vs without the floor) is
  still pending.

## LANDED — per-spectrum TSV + scoring progress in the eval scorer
`NuclideConfigEval`'s scorer is the `static score_config_over_precomputed` at `PeakFitImprove.cpp:1006`
(a *different* function from the same-named one in `NuclideConfig_GA.cpp`, which only runs during
`--eval-html` gallery render — after a ~1 hour 1464-chart render).  Added, gated by env
`PEAKFIT_PERSPEC_TSV=<path>`: one flushed row per spectrum (`source, detector, fg_cost, status,
miss_fraction`) plus a `[score] N of M` stderr counter.  The scoring phase previously had ZERO
progress visibility.  (An earlier copy mistakenly added to the `NuclideConfig_GA.cpp` gallery scorer
was reverted.)

## DIAGNOSIS — NaI low-resolution failure modes (default_config, 1464 spectra, R1.1 metric)
`total_cost=8025`, avg 5.48.  Cost concentrates in two modes covering ~76% of cost from ~16% of spectra:
- **Hard FAILED (fit ≠ Success): 49 spectra, cost 4900 = 61% of total.**  Clusters: **Pd103 ×11**,
  Th228 ×6, Lu176 ×5, Tl201wTl202 ×4, Sc46 ×4, Co58 ×4.
- **Finds-nothing (ok, miss_fraction ≥ 0.8): 187 spectra, cost 1167 = 15%.**  Clusters: Pd103 ×24,
  U235 ×16, Ta182 ×16, Pu239 ×13, Ra223 ×11, U238 ×10, Kr85 ×10, U233 ×9, Lu177 ×9, Ba133 ×9.
- Wins: 618 spectra at cost ≤ 0.
**Two-pronged interpretation:** (1) the 49 failures are genuine fit-robustness bugs worth debugging
(Pd103 low-E x-ray forest is the standout — 35 problematic Pd103 spectra total); (2) much of
finds-nothing is sources whose truth peaks NaI **physically cannot resolve** (x-ray forests, weak
single lines like Kr85's 514 keV @0.4% BR) — the metric penalizes the fitter for the unfindable,
which is the SAME truth-reliability issue as the scoring review (R1.2/R3).  So "better low-res score"
is partly a metric fix, not a fitter fix.  **Next:** debug the Pd103/Th228 failure cluster; fold
low-res truth-realism into the next scoring increment.

---

# 2026-07-17 — Phase 3 · R6: auto co-fit of strong unmodeled interfering lines

Fixes the standing [HIGH] algorithm item below ("A source whose line interferes with a strong
unmodeled NORM/background line must be co-fit").  Full plan:
`~/.claude/plans/you-are-continuing-phase-rippling-creek.md`.  Landing in three build-test-review
increments (user chose full R6 + eval/gallery accuracy pass after each behavioral increment).

## Increment 1 LANDED 2026-07-17 (working tree; NOT committed) — detection only, no fit-path change
- [x] **NORM table refactored to `{energy, parent}` pairs** (`src/FitPeaksForNuclides.cpp`): the flat
  `std::array<double,51> sk_strong_norm_gamma_energies_kev` became
  `std::array<StrongNormGammaLine,51> sk_strong_norm_gamma_lines` (`StrongNormGammaLine{double energy;
  const char *parent_symbol;}`), each line attributed to its `get_norm_sources` parent
  (Pb214/Bi214→Ra226, Th234/Pa234m→U238, Ac228/Tl208/Bi212→Th232, U235, Ra226, K40; the 7 Pb/Bi/Tl
  K-xrays → nullptr).  `is_near_strong_norm_gamma` now iterates `.energy` — **behavior byte-identical**
  (energy multiset + order verified equal to HEAD; the only behavior-sensitive consumer,
  `NuclideConfig_GA.cpp:417`, unchanged).  Parent symbols resolve to `Nuclide*` lazily (never at
  static init).
- [x] **Detection helper** `detail::find_strong_unmodeled_interferers(...)` (declared in the header
  `detail` namespace for unit testing; defined in a re-opened `detail` block after the table).  Returns
  `InterfererCandidate{energy, nuclide/*null⇒float*/, detection_z, from_background_search}`.  Foreground
  path only in inc1: for each in-range requested-source gamma within `sm_interferer_near_num_fwhm=2.0`
  ×FWHM of a strong NORM line whose parent is not already modeled and is not a requested source, and
  which is NOT explained by the source's own chain (`sm_interferer_source_owns_num_fwhm=1.0`), require a
  foreground auto-search peak within `sm_interferer_confirm_num_fwhm=0.5`×FWHM at
  area/uncert ≥ `sm_interferer_min_detect_z=5.0` (Poisson `sqrt(area)` fallback).  Doublet guard
  (single-line source vs single-line interferer <`sm_interferer_doublet_min_fwhm=1.0`×FWHM → skip+warn)
  runs BEFORE source-owns so it wins the warning.  Caller pre-expands `sources`→`RequestedSourceGammas`
  so the helper is DB-light and unit-testable on synthetic input.  Nullable `GlobalContinuumEstimate*`
  (forward-declared; R1-step2 stub — never dereferenced) reserved for the gross-counts confirmation.
- [x] **6 file-local `constexpr` `sm_interferer_*` thresholds** added to the safety-rail block (not
  GA-exposed for now; become the tunables to iterate in the accuracy pass).
- [x] **Unit test** `test_interferer_detection_unit` (StatisticalDetailHelpers suite): K40-for-Eu152
  positive; suppressed when K40 requested / unconfirmed / source-owns / fit_norm_peaks; doublet-guard
  warning; `is_near_strong_norm_gamma` probe-energy regression.  Full `test_fitPeaksForSources` green
  (24 cases; trinitite fits unchanged — no fit-path change yet).
- **Adversarial review (subagent, 2026-07-17): NO CONFIRMED BUGS.**  Independently verified: energy
  multiset+order identical to HEAD; all 44 parent attributions correct vs `data/sandia.decay.xml`
  (incl. the subtle U238-prompt→Th234/Pa234m vs Ra226-separate→Pb214/Bi214 chain split); helper margins
  load-bearing and correct (near 2.0 > owns 1.0 is why Eu152-1457/K40-1460 co-fits); guard ordering,
  boundary operators, null-checks, no UB.  Three low-severity carried notes:
  1. **Cross-source order dependency** (`emitted_line_energies` shared across sources): if two requested
     sources are both within 2 FWHM of the SAME single-line interferer (realistically only K40 1460),
     iteration order decides doublet-skip vs co-fit.  Rare; detection-only.  Left as a known limitation.
  2. **De-dup is by line energy, not nuclide** — a parent with two confirmed interfering lines yields
     two candidates.  **Increment 2 MUST de-dup by nuclide when building the interferer rel-eff curve**
     (add the nuclide once; widen a ROI per candidate line).
  3. **Partial-yields robustness**: FIXED inline — `num_significant_lines` now returns "multi-line"
     (conservative, no spurious doublet) when a caller passes `yields.size() != energies.size()`.
- **GA note:** none — no gene semantics changed; `is_near_strong_norm_gamma` byte-identical.  The
  `sm_interferer_*` thresholds are the manual-loop tunables for the accuracy pass after inc 2/3.

## Increment 2 LANDED 2026-07-18 (working tree; NOT committed) — nuclide-interferer co-fit wired in
- [x] **Solve integration** in `fit_peaks_for_nuclide_relactauto` (`src/FitPeaksForNuclides.cpp`): a
  `if( !fit_norm_peaks )` block after the NORM-curve block expands `sources`→`RequestedSourceGammas`
  (`get_source_photons`/`get_source_age`), calls `find_strong_unmodeled_interferers`, and for each
  nuclide candidate (de-duped BY NUCLIDE) builds ONE extra `RelEffCurveInput` cloned from the NORM
  recipe (FramPhysicalModel, order 0, soil self-atten; NORM parents get `get_norm_sources` equilibrium
  ages, non-NORM get age 0).  A covering ROI is ensured per interfering LINE, **merging any ROIs the
  coverage window touches into one** so `resolve_overlapping_rois` never sees a non-escape overlap
  (which it asserts against).  All best-effort under try/catch.  Tracking locals `auto_interferer_nucs`
  + `auto_interferer_float_energies` declared by the rel-eff-index trackers.
- [x] **Result filter**: a dedicated branch at the top of the result-populate loop drops peaks whose
  `parentNuclide() ∈ auto_interferer_nucs` (or source-less peaks matching a float energy), REGARDLESS
  of `norm_peaks_dont_use`, before the existing NORM branch; a `result.warnings` entry names what was
  co-fit.  Interferer curve + its ROI persist across no-ecal/desperation/refinement via the `options`
  copies and refinement re-clustering (verified by trace).
- [x] **Tests**: `test_trinitite_default_sequence` + `test_trinitite_do_not_use_existing_sequence`
  converted `{K40,Eu152}`→`{Eu152}` (+ assert no returned K40 peak); new
  `test_eu152_interferer_matches_joint` (Eu152-alone ≈ joint).  `test_eu154_does_not_remove_strong_eu152_peak`
  KEPT on joint `{K40,Eu152}` (regression: requested K40 disables the auto-interferer).  Full
  `test_fitPeaksForSources` GREEN (25 cases).
- **Empirical result (trinitite):** R6 fires (adds K40 + ROI `[1457.19,1464.45]`), but K40 fits to
  ~0 — the ~modest net 1460 peak is attributed to Eu152 1457 / continuum.  The **joint fit does the
  same (K40≈0)**, so Eu152-alone (R6) EXACTLY matches joint: `auto 16 Eu152/0 K40 == joint 16 Eu152/0
  K40`, key-line areas within tolerance.  R6 is correct + no regression here; trinitite is just not a
  strong *demonstration* (its 1460 line is weak/ambiguous).
- **Review status:** the increment-2 adversarial subagent hit an API session limit mid-run.  Its one
  flagged concern (post-filter bystander/observable logic re-introducing interferer peaks) was traced
  by hand: `compute_observable_peaks(result.fit_peaks)` starts from the FILTERED peaks and the
  bystander path only re-adds USER peaks (`updated_bystander` = copy of an existing user peak), so no
  interferer nuclide can be re-introduced (confirmed empirically: 0 K40 in observable).  ROI-merge,
  curve construction (FramPhysicalModel order 0; coexists with the non-physical sources curve exactly
  as the NORM path does), and the result-filter branch (never touches 511/escape floats) self-reviewed
  clean.  **Re-review the full diff when the session limit resets.**
- **OPEN (verify with evidence, not yet reproduced):** *observable-refit re-absorption for a STRONG,
  UNAMBIGUOUS interferer.*  Filter-early drops the interferer before `compute_observable_peaks`, whose
  per-ROI `fit_peaks_in_spectrum_LM` refit (and the earlier `combine_overlapping_peaks_in_rois`) could
  in principle re-attribute a strong interferer's counts to the source's nearby weak line.  Trinitite
  cannot exercise this (K40≈0 in the solve for both auto and joint).  Plan: hunt a strong-interferer
  spectrum in the accuracy pass; if re-absorption is confirmed, evaluate filter-late (keep the
  interferer through the observable refit, drop from the final vectors) — noting `combine` may itself
  need to be made interferer-aware.  GA note: no gene semantics changed.

## Increment 3 ATTEMPTED 2026-07-18 — ambient Cs137/Co60 sweep; DISABLED after a joint-fit regression
- [x] **Foreground detection refactored** into a shared `check_interferer_line(es, src_single_line,
  le, parent_symbol)` lambda so the loop CAN sweep multiple line tables.  Behavior-preserving for the
  NORM path (unit-test Cases A–F unchanged).  KEPT.
- [x] **Ambient line table** `sk_ambient_interferer_lines` ({661.657→Cs137, 1173.228/1332.492→Co60})
  + `interferer_is_single_line` counting ambient lines.  KEPT (harmless: the ambient count is 0 for
  the NORM parents, so NORM doublet-guard behavior is unchanged), but the ambient SWEEP call is
  **commented out**.
- [!] **REGRESSION FOUND → ambient sweep DISABLED.**  Enabling the Cs137/Co60 sweep collapsed the
  `{K40,Eu152}` JOINT fit to **0 returned peaks** (reproduced in ISOLATION, so not contention): the
  Eu152-ALONE fit stayed healthy (18 Eu152 peaks), but co-fitting an ambient Cs137 interferer
  alongside TWO requested sources destabilized the RelActCalcAuto solve so every peak became
  insignificant.  Disabling the sweep restores exact inc2 behavior
  (`test_eu152_interferer_matches_joint`: auto 16 == joint 16, GREEN).  Unit Case G removed with it.
  - **What this reveals (investigate before re-enabling ambient OR the deferred background/float
    paths):** adding an interferer rel-eff curve can DESTABILIZE a MULTI-SOURCE solve.  inc2's NORM
    co-fit passed all tests + the eval, and the single-source (Eu152-alone) ambient co-fit is fine, so
    the trigger is (extra requested source) × (extra interferer curve).  Hypotheses: the
    FramPhysicalModel order-0 interferer curve is under-determined from a single ambient line and,
    with a second source curve present, drives a degenerate solve; or the added Cs137 661 ROI
    starves/overlaps an Eu152 region.  **Fix ideas:** cap total interferer nuclides; use a simpler
    (LnY / constant) interferer curve for a single-line ambient; or add single-line ambients as a
    `FloatingPeak` instead of a nuclide curve.  **This is ALSO a latent risk for inc2** on untested
    multi-source spectra — a targeted multi-source stress test is warranted before committing.
- **Also deferred:** a dedicated background auto-search (disambiguation) + unattributable FLOATING
  interferer peaks (5-site re-apply helper).  The helper still accepts `background`/`drf`/
  `peak_fit_prefs` (reserved).
- Full `test_fitPeaksForSources` GREEN with the sweep disabled (= inc2 + a cosmetic refactor).

## Accuracy pass (2026-07-18) — R6 isolated, HPGe, 127-spectrum holdout sample, default_config
Method: `PeakFitImprove --action NuclideConfigEval --det-type HPGe --config-genes default
--live-time 300_seconds --holdout-frac 0.03 --holdout-role eval` (avg cost; more-negative = better).
**A confound was found and controlled:** the first "baseline" binary predated the committed
extent(eba205a4)/SNIP(45051549) work, so a clean R6-off was rebuilt from HEAD via a scoped
`git stash` of ONLY the two R6 files (restored immediately; verified my +483/-31 diff intact).
- **R6-off (HEAD, extent+SNIP, no R6): avg cost −3.656**
- **R6-on (HEAD + R6 inc2): avg cost −5.101**  ⇒ **~40% better**, isolated to R6.
- The clean R6-off exactly reproduced the pre-extent baseline (−3.656), confirming determinism and
  that extent/SNIP don't move this metric.  R6 improves source-area accuracy by co-fitting NORM
  contamination near source lines (co-fit peaks are dropped from results, so no extra peaks are
  reported — the gain is cleaner source areas).  **Caveats:** small sample (127 HPGe spectra, 3%
  holdout), UN-tuned `default_config`, aggregate metric; the SpectrumChartD3 galleries
  (`run_r6_eval/baseline_hpge_r6off.html`, `r6on_inc2_hpge.html`) and a full-dataset eval (big
  machine) + per-class `sm_interferer_*` tuning remain for user review.  Machine here is 10-core, so
  the full eval was not attempted.

---

# 2026-07-17 — Phase 3 (architectural).  R1 step 1: SNIP continuum rewrite + FWHM-window estimator

Phase 3 restructures the gate-then-fit pipeline (see the session plan doc
`~/.claude/plans/you-are-working-on-tingly-abelson.md`).  Ordering per user: SNIP continuum FIRST
(foundational + a licensing win), then manual-tuned baseline, R6 interferer co-fit, R1-step2 global
continuum threading, R2 rescue.

## Landed 2026-07-17 (working tree; NOT committed)

- [x] **SNIP clean-room rewrite + ROOT/TSpectrum removal** (`src/PeakFit.cpp`, `InterSpec/PeakFit.h`,
  `NOTICE.html`).  The ROOT-derived `calculateContinuum` (~770 lines) and dead
  `findPeaksByRelaxation` (~450 lines, no callers) were DELETED, along with the `kBackOrder*/
  kBackSmoothing*` enum.  This removes the last CERN-ROOT-derived (LGPL) code from the project —
  a step toward the eventual MIT relicense (Minuit2, still LGPL, is a separate remaining item).
  NOTICE.html's ROOT/TSpectrum credit removed (Minuit2 entry kept).
  - New file-static `calculate_continuum_variable_window(spectrum, ssize, channel_windows,
    filter_order, presmooth, lls)` implements SNIP from the published algorithm (Ryan et al. NIM B
    34 (1988); Morhac et al. NIM A 401 (1997)): vanishing-central-difference clip estimates
    E2/E4/E6 (binomial coeffs), max-composition, `min(prev, est)`; per-channel increasing-window
    schedule `d_j(m)=min(m, w_j)`; optional 3-ch boxcar presmooth; optional LLS transform.
  - Legacy `estimateContinuum(data)` reimplemented on the new kernel with a constant 125-channel
    window, order 6.  **Numerically validated equal to the deleted ROOT code to max rel diff
    ≤2.5e-6 over 8 detector configs × 223 spectra** before deletion.  Tiny-spectra path (<5 ch)
    now returns data-as-continuum without asserting (old code silently no-op'd there).
- [x] **FWHM-window estimator overload** `estimateContinuum(data, fwhm_at_energy, num_fwhm_window,
  filter_order, presmooth, lls)`: per-channel clip half-width `w_j = round(k·FWHM(E_j)/
  channel_width_j)`, clamped `[2, nchan/2-1]`, nearest-valid-neighbor fill.  Addresses the user's
  finding that the fixed 125-ch window (~60 keV on HPGe) is far too wide.
  - **Key finding (quantified):** at FWHM-scaled windows the order-6 composite's E4/E6 taps
    (w/2, w/3, 2w/3 from center) land INSIDE a peak of that width and pull the estimate up,
    leaving a residual under strong peaks — over the 223-spectrum HPGe-30s set, order-6 at
    k=1.5 preserves only `frac(z>3)=0.001` of peak bins vs `0.059` for order-2 (and `0.059`
    for the legacy wide window).  **Order 2 is the right filter for FWHM-scaled windows; order 6
    only pays off when the window ≫ peak width (the legacy 125-ch case).**  This is why the
    overload defaults callers to order 2.
  - Sweep assessment (median (data−snip)/√data and peak-preservation `frac z>3`, per energy
    third, all 8 configs): order-2 k≈1.0–1.5 tracks the continuum in flat regions (median z≈0)
    while preserving peaks like the legacy window, at a physically-scaled window ~10× narrower on
    HPGe.  Wider k and LLS pass further under peaks but dip lower in structured low-E regions;
    presmooth slightly reduces peak-standout.  **Pending user sign-off on the refreshed viewers.**
  - Tooling: `target/peak_fit_improve/SnipContinuumPlots.cpp` extended (per-config generic FWHM
    fns; 8 configs = {Detective-X HPGe, R500-NaI, IdentiFINDER-LaBr3, CZT H3D-M400} × {30,300 s};
    7 swept variants k1.0/1.5/2.0/3.0-o2, k1.5-o6, k1.5+smooth, k1.5+LLS).
    `snip_continuum_plots/make_viewer.py` rewritten multi-series (toggleable legend).
  - Test: `test_estimate_continuum_snip` in the StatisticalDetailHelpers suite (min-filter
    property; order-2 peak erasure; order-6 under-peak residual; legacy-equivalence via a
    constant-125 window at order 6; LLS/presmooth finiteness; invalid-input throws).
  - **Inline adversarial review** (the subagent hit a session limit; done by hand): E2/E4/E6
    coeffs verified against the 2p-th central-difference derivation; `min(prev,est)` proven never
    to exceed data (cross-iteration, since prev is monotone non-increasing); all E4/E6 taps proven
    within `[j-w, j+w]` (3·⌊w/3⌋≤w, 2·⌊w/2⌋≤w) so no OOB; overload window clamp proven to satisfy
    `2·w_max+1 ≤ nchan`; LLS inverse verified exact (float rel err ~2e-4 only at y~1e8).  One bug
    found and fixed: the legacy wrapper asserted/failed on <3-channel spectra (window forced to 1
    exceeded `2w+1≤nchan`); now guarded.
  - **GA note:** none yet — no production consumer of the new estimator until R1-step2 threads it
    into the gates.  The `filter_order`/`num_fwhm_window` choices become tunable inputs there.
- [x] **Extent-restricted continuum + in-app ROOT credit removal** (2026-07-17, follow-up to user
  review of the viewers).  User feedback: order-2 k=1.0 best for low/med, k=1.5+smooth best for
  HPGe, both weak below ~230 keV where the detector turn-on dominates; wants ONE variant for all
  classes.  Two changes:
  - `estimateContinuum(...)` FWHM overload gained optional `restrict_lower_energy/
    restrict_upper_energy`; the kernel gained a `[clip_begin, clip_end)` sub-range.  Channels
    outside the range are left equal to the data (zero gating excess), and taps never read across
    the boundary — so the sub-extent turn-on hump cannot pull the in-range continuum up.  Callers
    pass the valid spectroscopic extent from `FitPeaksForNuclides::find_valid_energy_range()` (the
    resolution-aware `LowerExtentFinder`, committed this same day in `d8dd92db` by the concurrent
    extent session — this SNIP work sits directly on top of it).
  - **Unified candidate = `k1.5sx`** (1.5×FWHM, order 2, presmooth, restricted to extent): the
    HPGe favorite plus the extent fix that recovers the low/med sub-230 region.  Default-on in the
    viewers; awaiting user confirmation it serves all classes.
  - In-app ROOT/TSpectrum credit removed from `InterSpec_resources/static_text/
    copyright_and_about.xml` (About → Credits tab), matching the NOTICE.html removal.
  - Test extended: `test_estimate_continuum_snip` now asserts out-of-range channels equal the data
    and the in-range peak is still clipped.  Full StatisticalDetailHelpers suite green.
  - **Zoomable review pages** (`SnipContinuumPlots.cpp` `write_snip_d3_html`): one self-contained
    `snip_d3_<config>.html` per detector×live-time, built with the real `D3SpectrumExport` /
    SpectrumChartD3.js (mirrors `FitPeaksForNuclideDev::write_spectroscopic_extent_html`) — each
    spectrum is a zoomable InterSpec chart with the foreground + the SNIP continuum (as a
    Background `Measurement`) + orange extent reference lines.  D3 pages inline all JS (portable
    but large; HPGe ~85 MB/page).  Also emits the lightweight custom-canvas viewers
    (`snip_<config>.html`) and `snip_gallery.html`.
- [x] **SNIP tuning resolved to TWO per-class parameter sets, keyed on the existing HPGe-vs-NotHPGe
  classification** (2026-07-17, round-2 user review of the D3 pages).  User flagged two residual
  artifacts on the round-1 `k1.5sx`: HPGe continuum rides up into tall peaks (cosmetic); CZT
  continuum undershoots at high energy.  Quantified fixes over the 223-spectrum sets:
  - `presmooth` generalized from a bool to an int half-width (0/1/2/3 = off/3ch/5ch/7ch boxcar) in
    the kernel + overload; the legacy wrapper passes 0.
  - **Finding: no single (window, smooth, LLS) works for all classes.**  LLS fixes the HPGe
    tall-peak pedestal (Ag110m 884 keV under-peak ratio 2.89→1.24 at k2.0+LLS) but *wrecks* NaI
    Compton-edge tracking (Cs137 400-600 keV q75 c/data 0.98→0.40, because LLS compresses the
    Compton step so the min-filter clips across it).  Extra smoothing fixes the CZT sparse-high-E
    undershoot (At211 300-650 keV q75 0.50→1.00 at 7ch) and does NOT regress NaI/LaBr (0.98→1.00),
    but makes HPGe worse without LLS.  The split falls exactly on HPGe vs not-HPGe:
    - **HPGe (High):** `2.0×FWHM, order-2, 3ch presmooth, LLS on, restricted to extent`.
    - **NotHPGe (NaI/LaBr/CZT = LowOrMedRes):** `1.5×FWHM, order-2, 7ch presmooth, LLS off,
      restricted`.  One set serves all three low/med classes (verified NaI, LaBr, CZT).
  - This is consistent with the rest of the pipeline: the estimator is ONE function whose
    (window, smooth, LLS) are selected by `PeakFitUtils::coarse_det_type` — the same 2-way
    HPGe/NotHPGe split (96.6% accurate) that already picks the FWHM curve.  These become the
    GA-tunable inputs in R1-step2 (start from these two sets).  D3 pages now draw each config's
    class-appropriate variant (HPGe→k2.0slx, low/med→k1.5s7x); the sweep JSONs carry all
    candidates for comparison.  **Pending user sign-off on the two per-class continua.**
  - **Concurrency note:** the extent session committed `d8dd92db` (415-line change to the SAME
    `src/PeakFit.cpp` I am editing) into the shared working tree mid-session; verified my
    uncommitted SNIP diff and their committed `LowerExtentFinder` coexist with no reversion.  Their
    activity (a `git clean`, most likely) deleted the untracked `snip_continuum_plots/` output dir
    — `make_viewer.py` was recreated from context; the README and hand-made reference PNGs
    (`s_cs137.png`, `s_nai_th232.png`) were untracked and are lost (regenerable output only).

---

# 2026-07-12 review round (post genome-v2)

Second multi-agent review focused on (a) the genome-v2 machinery (adaptive sideband extent,
clean-gap merge, AICc continuum selection — added 2026-07, not covered by the 2026-06 round),
(b) fresh bugs in the pipeline, (c) approach-level critique.  High-severity claims were re-verified
by direct trace (two initially-reported items were REFUTED on re-trace and are recorded below so
they aren't re-raised).  Full write-up with statistical analysis and Phase-2/3 improvement plan:
see the session plan doc (findings summarized here in punch-list form).

## Fixed this round (2026-07-12/13)
- [x] **(P1) FitNormBkgrndPeaks path discarded existing-ROI setup** — `fit_peaks_for_nuclide_relactauto`'s
  NORM branch rebuilt `options.rois` from raw `input_rois`, throwing away the existing-ROI trims and
  mixed-ROI fixed-bounds ranges (bystander floats then pruned by `remove_floating_peaks_without_roi`).
  Now rebuilds from the current (trimmed) `options.rois`.  Also deleted a DEAD full manual rel-eff
  solve (`norm_rois`, never read; the caller already merged NORM ROIs) — wasted CPU on every
  NORM-enabled GA eval — plus its now-orphaned `norm_drf`/`norm_config`/`ext_atten` setup.  Note two
  intentional side effects of the deletion: (a) if that pre-solve threw, the old code silently
  skipped adding the NORM rel-eff curve for the whole fit — the curve is now always added; (b) the
  concrete external attenuator only ever fed the dead pre-solve and was never applied to the NORM
  curve — a comment now records where to add it if wanted.
- [x] **(P2) compute_filtered_chi2_dof mixed calibration frames** — it paired TRUE-energy
  `m_final_roi_ranges` with SPECTRUM-cal `m_peaks_without_back_sub` against the caller-supplied
  (possibly cal-advanced) loop `foreground`; whichever spectrum was passed, one side was displaced by
  the fitted cal shift (worst: in-loop accept gate scored the incumbent with peaks displaced by the
  full just-applied shift while the challenger was scored nearly consistent; final filter wrong by
  full shift on rejection/failure exits and under DoNotRefineEnergyCal).  Now renamed
  `compute_filtered_chi2_per_channel` (it always was chi2/channel, not chi2/dof) and internally
  consistent: `m_final_roi_ranges_in_spectrum_cal` + `m_peaks_without_back_sub` +
  `solution.m_foreground`.
- [x] **(P3) Desperation retry lost better solutions / re-enabled ecal** — the Physical-Model
  desperation solve overwrote the no-ecal `trial_solution` (a Success-but->10 no-ecal fit was lost if
  desperation failed or was worse), and `desperation_opts = options` silently kept
  `energy_cal_type = NonLinearFit`.  Now: no-ecal and desperation kept separately, final selection is
  best-of {ecal, no-ecal, desperation} by reduced chi2, and desperation runs with NoFit ecal.
  (Still gated on ecal-active per the outer condition; un-gating is a design question, deferred.)
- [x] **(P5) InitialRoi::estimated_amplitude had two meanings** — `estimate_initial_rois_without_peaks`
  stored a dimensionless br*eff yield where the other callers store peak-area COUNTS, silently
  disabling `find_clean_gap_between`'s tail test in the fallback path.  Unit now documented as counts;
  the no-activity fallback passes 0 = explicit "unknown" (documented gap-width-only behavior).
- [x] **(P6) Sibling-expansion resurrected dropped peaks** — the default-mode loop that removes user
  peaks sharing a continuum with a removed peak re-added a stale copy for EVERY sibling, including
  same-source peaks the fit had evaluated and found unobservable.  Now: same-source (in-fit) siblings
  are removed without a replacement copy; only other-source bystander siblings get the copy.
- [x] **(P7) FWHM-fit exception discarded a valid caller DRF** — with >6 auto peaks the peak-based
  FWHM fit is preferred even over a valid DRF; on exception the fallback replaced `drf` with a
  GENERIC detector (also corrupting downstream efficiency use).  Now falls back to the input DRF's
  own resolution info when present, generic only otherwise.
- [x] **(P8a) low_energy_analysis_floor comment said 30 keV non-HPGe; code uses 25** — comment fixed to 25.
- [x] **(P8b) Under-min-width ROI at emission silently dropped keep-gate-passed gammas** — in
  `cluster_gammas_to_rois` ROI emission, a cluster below `min_fwhm_roi` (via the `previous_roi_upper`
  clamp or a clean-gap/valley split sliver) was `continue`d away with its gammas.  Now folded back
  into the previous ROI when abutting (mirrors merge_rois' under-width guards); still dropped when
  isolated (essentially unreachable: the adaptive-extent core alone is ~2*roi_core_num_fwhm wide).
- [x] **(P8c) estimate_local_continuum silently relocated sidebands at spectrum edges** —
  `eqn_from_offsets` clamps sidebands inward near channel 0 / last channel, anchoring the line at the
  wrong place while `valid` stayed true.  Now returns invalid when a sideband would extend past
  either spectrum edge.
- [x] **(P8d) GUI removePeaks could throw uncaught on async model drift** — PeakModel::removePeaks
  throws on the first pointer-identity miss; the fit is async so the user can edit/delete peaks
  before Accept.  All three accept sites now use a tolerant per-peak remove
  (`remove_peaks_tolerantly`).
- [x] **(P8e) GUI preview showed old + replacement peaks stacked** — preview never subtracted
  `original_peaks_to_remove`; now previews (current − to-remove) + observable.
- [x] **(P8f) FitNormBkgrndPeaksDontUse + empty sources returned Success with zero peaks silently** —
  now pushes a warning when all fit peaks were NORM-filtered.
- [x] **(P8g) Bystander FloatingPeakResult matching could double-bind** — nearest-energy match with a
  0.5 keV tolerance could bind two close bystanders to the same result; each result is now consumed
  at most once.
- [x] **(P8h) compute_roi_chi2_significance p-value unguarded at 1.0** — when a ROI's chi2 improvement
  is ~0, `cdf(complement(chi2_dist, x))` rounds to exactly 1.0 and `boost::math::quantile(normal, 1.0)`
  throws an erfc_inv Overflow that propagated up and FAILED the whole fit (surfaced by the P2
  frame-consistency change perturbing a ROI into this corner in test_am241_smoke; latent before).
  Now clamped: p > 1-1e-12 → equivalent_z = -40 (insignificant), mirroring the existing p < 1e-300
  → +40 guard.
- [x] **(P4 — RETRACTED, no bug)** `DoNotRefineEnergyCal` does NOT return peaks in the fitted cal:
  returned peaks come from `m_peaks_without_back_sub` (spectrum cal), and the working foreground is
  never cal-advanced in this mode, so spectrum cal == original cal and the PeakFitResult contract
  holds.  Header doc note added.  (The mode's final-filter window mismatch was the P2 fix.)
- [x] **(RETRACTED, no bug)** merge_rois mid-FWHM interpolation "extrapolates / trips assert":
  with ROIs sorted by lower energy and overlapping, `mid_dist` is provably in [0,1] (nested case
  yields exactly 0); the assert cannot fire.  Optional clarity nit only (direct
  peakResolutionFWHM(mid_energy) is stubbed in comments).

## Phase 2 — statistical hardening: IMPLEMENTED 2026-07-13 (GA RE-TUNE NOW REQUIRED)
All items below changed tuned-threshold semantics; the per-det-type NuclideConfig GA must be
re-run before trusting any tuned config (they were already invalid per the objective-sign fix).
- [x] **estimate_local_continuum upgraded** (the shared fix for several items): optional
  caller-supplied `predicted_signal(x0,x1)` subtracted from each sideband sum; sideband windows
  slide outward (up to 3 widths) past interfering unfit auto-search peaks, else invalid; sideband
  counts/extents exposed on `LocalContinuumEstimate` + new `sideband_asymmetry_z()` helper.
- [x] **extend_roi_by_sidebands de-biased**: per-block threshold now Bonferroni-split so
  `extend_z` is the FAMILY-wise z for a full side (fixed the ~28% false-stop on flat continua at
  z=2 with ~7 blocks/side); block-z denominator now includes the estimation variance of the
  extrapolated anchor line (leverage-aware, from raw anchor counts).  Anchor-correlation of
  successive blocks noted as approximate (GA absorbs residual).
- [x] **find_clean_gap_between de-biased**: contamination test moved from per-~0.25-FWHM-block to
  WINDOW level (per-block understated window contamination by ~sqrt(block/window) ⇒ split bias);
  continuum anchor now signal-subtracted via the upgraded estimator (strong close peaks no longer
  inflate c_est and spuriously pass).
- [x] **merge_rois anchors grounded**: `InitialRoi::center_energy` now uses the largest
  auto-search peak's MEAN when one exists in the ROI (was always the geometric ROI midpoint), so
  the clean-gap groups and split constraints model signal at a real peak position.
- [x] **Keep-gate window = fit window** (also resolves #203): S/√(S+B) now evaluated with B over
  the cluster's CORE extent (outermost gammas ± roi_core_num_fwhm·FWHM — the always-included part
  of the emitted ROI), signal-subtracted sidebands, unfit-peak-relocating.
- [x] **Step continuum: probe geometry replaced by the agreed gated "let the model choose" ladder**
  (design recorded below): Gate 1 = z ≥ step_cont_min_peak_significance dominance (unchanged gene,
  B now signal-subtracted); Gate 2 = sideband_asymmetry_z ≥ 1.0 (sm_step_trial_min_asym_z, loose,
  cannot self-veto); Stage 3 = `trial_step_continuum`: fit_amp_and_offset (linear LS) trial of the
  selected polynomial vs its EQUAL-parameter-count step variant (Linear↔FlatStep,
  Quadratic↔LinearStep, so AICc penalties cancel), step kept only if it wins chi2 by
  `step_trial_chi2_margin` (new GA gene REPLACING step_cont_left_right_nsigma; genome size
  unchanged; NuclideConfig_GA serialize/parse/init/mutate/crossover updated.  NOTE: checkpoints
  predating the rename do NOT resume — from_string fails on the missing new key, the individuals
  are discarded, and an all-old checkpoint exits with an error.  Fail-loud is intentional: a
  population tuned against different gene semantics should not silently seed a new run).
  Trial peaks = cluster gammas merged within 1σ,
  capped at 8 largest (line-dense-source conditioning, per user feedback: amplitude-aware
  subtraction, no position-based neighbor vetoes).  Trials use NoSkew for both candidates (fitted
  skew params unknown at clustering time; identical shape keeps the comparison fair).
  `should_use_step_continuum` deleted.  Also fixes the old LinearStep/FlatStep ≥4-FWHM wart and
  the step-gate S/B window mismatch (B now signal-subtracted, dominant-window based).
  NOTE: Gate 2 requires a VALID sideband estimate — a ROI at a spectrum edge (or with hopelessly
  contaminated sidebands) never gets a step; previously the invalid-fallback path could still
  reach the probe test.

  **Known/intended behavior deltas of the Phase-2 statistics** (from the adversarial verify pass —
  interpret re-tuned thresholds against these semantics):
  - Window-level clean-gap testing relaxes detection of LOCALIZED in-gap contamination by up to
    ~sqrt(window/block) ≈ 2x in z (only the diffuse-tail direction was the motivating fix): a
    modeled in-gap line with amplitude between z·sqrt(B_block) and z·sqrt(B_window) that the old
    per-block test rejected now counts as clean.  Correct for the test's purpose (bias of a
    whole-window continuum anchor scales with total contamination), but a real relaxation.
  - The gross-count fallbacks (clean-gap when the continuum estimate is invalid; keep-gate
    b_est=data_area) still re-include the signal in the continuum estimate — now only reachable at
    spectrum edges / fully-contaminated sidebands, and conservative in direction.
  - Keep-gate "window = fit window" is approximate: the core uses the seed gamma's FWHM at both
    ends and omits the skew low-side allowance the fit extent adds; S covers >=92% of the core's
    Gaussian mass at the GA-minimum core width (few-% z overstatement).
  - `InitialRoi::fwhm` is still evaluated at the geometric ROI midpoint even when center_energy
    was re-anchored on a real peak mean.
  - The Bonferroni per-block threshold has a floor (~1.4-1.5 for ~7 blocks, where family alpha
    saturates at 1) and a discontinuity at extend_z==0; the GA gene range [0.75, 4.0] keeps both
    unreachable in practice.
  - merge_rois invalid-after-split handling changed (post-review fix): when the split leaves the
    current ROI invalid, it now folds back into the previous ROI when the combined width is
    acceptable (coverage preserved; the clean gap goes unused), or restores the previous ROI's
    original upper bound and drops current when merging would exceed the max width (previously:
    always dropped current AND left the previous ROI shrunk to the split point).

## Open — smaller items from this round
- **Mixed-ROI fixed bounds not preserved across refinement iterations** — header contract
  (FitSrcPeaksOptions doc: mixed ROI "will maintain the original bounds") is violated after the first
  accepted refinement re-cluster; bystander float can lose its covering ROI.  Fix: re-append mixed
  fixed ranges each refinement iteration.
- **rois_are_similar convergence test** compares refined ROIs vs the solution's INPUT rois, not
  `m_final_roi_ranges` (what the solve actually used).
- **Keep-gate erases cluster members before the keep decision** (`cluster_gammas_to_rois`): gammas
  absorbed into a REJECTED cluster can never seed their own ROI.  Possibly intentional (satellites
  are weaker than the seed) — needs a comment or a fix.
- **Existing-ROI trimming duplicated** with divergent rules (pre-fit trimmer vs
  filter_rois_for_existing: different abut rules / overlap buffers); NORM nuclide list hard-coded at
  5 sites; generic-DRF switch at 5 sites.  Consolidate.
- **Three lossy output vectors** (`uncombined_fit_peaks` → combine → `fit_peaks` → filter/refit →
  `observable_peaks`): source attribution lost before the only vector GUI/GA consume; bystander/
  sibling reconstruction bolts un-combined peaks onto the combined vector.  Design fix: one peak list
  with provenance/observability flags (Phase 3).
- **Observability via independent LM refit** whose continuum/bounds can diverge from the RelActAuto
  solve; significance computed 3 different ways; ROI bound-expansion hack papers over mean migration.
- **Manual-stage pruning is irreversible** (one-way information flow): a line dropped as a
  contaminant / by the keep-gate can only return via the narrow edge-ROI-restore heuristic; the auto
  stage never re-tests dropped regions against its improved rel-eff (Phase 3: fit-then-prune).

---

## Already fixed (this work, 2026-06-03/04)
- [x] **GA objective sign inverted** — reward terms summed into a minimized cost (PeakFitImprove.cpp ga_eval); negated.
- [x] **GA failure penalty 100.0 off-scale** — now documented `sm_fit_failure_penalty`.
- [x] **GA seeds from default_config(Low) for all det types** — now seeds from `sm_base_det_type`.
- [x] **Forced `local_debug_printout=true` (data race / I/O flood)** — removed.
- [x] **Init rel-eff eqn_order ignores #nuclides** — added dof-aware order clamp.
- [x] **Fallback activity live_time inconsistency** (estimate_initial_rois_fallback) — removed spurious `*live_time`.
- [x] **Desperation chi2 gate tested `solution` not `trial_solution`** + dof==0 divisions — added `reduced_chi2()` helper, fixed operand.
- [x] **FWHM-fit catch set got_fwhm_fcn=true** (skipped generic fallback) — now false.
- [x] **Cal-consistency**: translate/observable used loop-local `foreground` instead of `solution.m_foreground` — now `solution_foreground`.
- [x] **(new) Manual rel-eff acceptance metric** — was stat-only chi2 (always fell back for strong sources); now uses a systematic floor + excludes peaks the source can't account for, via stored `RelEffSolution::m_predicted_peak_counts`.
- [x] **(new) compute_observable_peaks refit kept inflated false-positives** — the per-ROI refit used `refitPeaksThatShareROI_LM`, which returns an EMPTY vector when a peak becomes insignificant; the caller could not distinguish that from a numerical failure and so KEPT the original (inflated continuum-curvature/degenerate "peak", e.g. the no-bg Eu152 1347/1363 false positives). Now uses `fit_peaks_in_spectrum_LM`, which reports insignificant peaks in `lost_peaks` so they are DROPPED; only a genuine `status != Success` falls back to keeping the current peaks. ("Option A".)
- [x] **(new) K40 must be fit with Eu-152 on trinitite** — Eu-152's 1457 keV line sits on K40's strong 1460 keV line; fitting Eu-152 alone mis-attributes that ROI and skews the rel-eff curve. Trinitite unit tests now fit `{K40, Eu152}` jointly (test_trinitite_default_sequence, test_trinitite_do_not_use_existing_sequence, test_eu154_does_not_remove_strong_eu152_peak). See open algorithm item below for the general (non-test) fix.
- [x] **(#88) FWHM-trim clamps unguarded** — the four `std::clamp(<energy>, fwhm_lower_energy, fwhm_upper_energy)` sites in `fit_peaks_for_nuclide_relactauto` (overlap-trim + two existing-peak-margin loops) now carry the same `have_fwhm_range = (lo>0 && hi>0 && lo<hi)` guard used everywhere else, so a `-1`/degenerate `[lo,hi]` no longer evaluates FWHM at a bogus energy or risks inverted-range UB.
- [x] **(#118) Rel-eff extrapolation — UPPER side now DRF-shaped** — both the auto-stage clamp (`fit_peaks_for_nuclide_relactauto`) and the manual-stage `manual_rel_eff` keep a flat clamp on the LOWER side but, above the fitted/matched span, shape the boundary rel-eff by the DRF intrinsic-efficiency falloff (`rel_eff(E)=rel_eff(E_hi)*drf_eff(E)/drf_eff(E_hi)`). New anon-namespace helpers `generic_drf_for_rel_eff_extrap` (null→generic by det_type) and `shape_rel_eff_above_boundary` (self-clamping via `intrinsicEfficiency`, degrades to flat hold when DRF unusable). Restores realistic high-E droop instead of an over-optimistic flat hold.
- [x] **(#103) compute_observable_peaks cal mismatch — VERIFIED RESOLVED** (no code change). The `solution_foreground` fix made it self-consistent: `m_foreground` stays in original/spectrum cal, `result.fit_peaks` (from `m_peaks_without_back_sub`) are also original cal, `translate_peaks_to_orig_cal` runs AFTER the call, so the significance window and `fit_peaks_in_spectrum_LM` refit both see peaks and data in the same cal. (Minor: stale comment "refit on fitted-cal foreground" at the call site is misleading — `solution_foreground` is original cal.)
- [x] **(#113) Source-blind peak merging — ASSESSED BENIGN** (no change). `combine_overlapping_peaks_in_rois`→`combine_peaks` keeps only the dominant source, and the merged set does flow into `observable_peaks`; BUT all GA scorers match by energy/ROI only and never read `parentNuclide()`, so the dropped attribution can't change a score. The only real (pre-existing, scorer-level) effect is two close distinct lines counted as one — a property of the energy-window matching, not of the merge. Leave as-is (optional: note it could matter for any external consumer that relies on `parentNuclide()` of returned peaks).
- [x] **(#158) Misleading 'best RelEff curve type' messaging** — dev-driver log strings + header doc now say a single configured curve type is fit (with an optional Physical-Model desperation retry), instead of implying a multi-curve ranking that does not exist.
- [x] **(#178a) Insignificant-ROI filter cal mismatch** — `insignificant_roi_ranges` now built from `m_final_roi_ranges_in_spectrum_cal` (was true-energy `m_final_roi_ranges`), so the 1 keV bound match is spectrum-cal vs spectrum-cal; only changes NaI/CZT cases where energy-cal is fit and shifts >1 keV. (Part (b), significance cal-consistency, deferred — needs GA awareness.)
- [x] **(#193) Zero amplitude-uncert hard-fails the manual solve** — manual seed now substitutes a Poisson `sqrt(max(area,1))` when `amplitudeUncert()<=0`, so a peak without a fitted uncert no longer trips RelActCalcManual's zero-uncert throw. (Same `-1` sentinel class as #238.)
- [x] **(#253) FramPhysicalModel order>0 debug-assert** — `genes_to_settings` now forces `rel_eff_eqn_order=0` for FramPhysicalModel, so the GA never samples the combination that aborts a debug build (fitness-equivalent in release).
- [x] **(#258) Non-deterministic best-individual fg/bg re-eval** — `SO_report_generation` now reads the cached `middle_costs.objective_fg/bg` of `best_chromosome_index` instead of re-evaluating, removing the extra full-dataset eval and reconciling exactly with `best_total_cost`.
- [ ] **Deviation-knot cap** — explored, then REVERTED (not a reliable speed lever; Ceres-convergence-dominated).
- All `test_fitPeaksForSources` cases pass with the above (Trinitite suite re-run green after the #88/#118/#178a/#193 changes).
- [x] **(#238) `-1.0` uncertainty sentinel mishandled** — added INLINE guards (no helper) at each significance/uncert site: escape S.E./D.E. significance (`assign_escape_peak_relationships`), the final-significance drop in `compute_observable_peaks`, and `combine_peaks` (guarded `eff_uncert` before squaring so `-1` no longer becomes `+1`). Each falls back to Poisson `sqrt(area)` when `peakAreaUncert()<=0` and `area>0`, else significance 0. (The cross-subsystem feed at 5416 was already handled by #193.) Left debug/serialization prints honest about "unset". PeakFit/PeakFitLM `999.0`/`0.0`/unguarded divisions left as-is (out of #238 scope).
- [x] **(#198) fallback `merge_rois` got null foreground** — `estimate_initial_rois_without_peaks` now takes and forwards `foreground` (so the no-matched-peaks path uses real spectrum-valley splits + channel-width gaps instead of midpoints). Intentionally still passes `{}` for unfit peaks (this branch is entered because nothing matched, so forwarding auto peaks as "unfit" would over-suppress merges). NOTE: `foreground` is effectively never genuinely null at the top level — the null was purely the un-forwarded default.

## Open LOW items (triaged 2026-06-06) — see deep-dives in session notes
- **#163** FitNorm nulls a supplied background (decision: safe doc+warning vs change `…DontUse`). **#168** fwhm_estimation_method not GA-tunable (Phase-4). **#183** edge-ROI restore heuristics. **#188** ROI significance veto Quadratic-only null (GA re-tune). ~~**#203** cluster keep-window vs emitted-ROI mismatch~~ (RESOLVED by the 2026-07-13 Phase-2 keep-gate core-window change). **#208** source+NORM gammas pooled before clustering. **#223** compute_significance_in_region S/√(S+B) saturation. **#248** interferer exemption uses a global gamma list (also the merge-prevention twin ~3958-4026); fix = scope to overlapping-ROI gammas (GA re-tune).
- **#173 (RECONSIDERED — do NOT count ROIs):** at the AUTO stage the bound is NOT ROI/line count — a single WIDE ROI (e.g. Pu ~600 keV) holds several gammas at distinct energies and can both support and need a curved multi-coefficient rel-eff, so a count clamp would wrongly flatten it. The only hard auto-stage bound is rank-based (can't fit N+1 coeffs from < N+1 distinct line ENERGIES); statistics/conditioning beyond that is the GA's job. Decision: leave auto order to the GA (+ RelActCalcAuto's internal manual-seed reduction); if ever needed, bound by distinct line energies, not ROIs. The MANUAL stage (count-peaks-vs-nuclides, "pretty close") is already guarded by the landed dof clamp. → effectively GA-controlled / defer.
- **Skip/defer:** #213 (weighted-mean re-derivation already implemented), #218 (minor √gross+floor nit), #228 (NORM⇄subtract exclusivity, design), #243 (combine_peaks width/uncert is reporting-only/benign).

## Known limitations / decisions
- **GA RNG is wall-clock seeded** (openGA `rng` is private, no setter) → runs not bit-reproducible. Per-genome scoring IS deterministic. Needs vendoring/patching openGA to fix (deferred).
- **No-background (config A) fits are slow** on NORM-/low-energy-dense sources (Ceres convergence), e.g. am241_smoke ~314 s. Deprioritized.
- **~~High-E weak-line loss~~ (RESOLVED — was the opposite)**: the Eu152 1347/1362 "lines" are FALSE POSITIVES — there are genuinely no peaks there (user inspected the spectrum). They were being *retained* (not lost) in the no-background path because the final per-ROI refit returned the input peak when it became insignificant (now fixed via "Option A", above). The bg-subtracted path correctly omitted them. The real underlying cause was the K40/Eu152 1457↔1460 interference (see algorithm item below).

---

## All confirmed findings (by severity)

NOTE: this is the full confirmed list; items resolved this session are checked off in
"Already fixed" above and also still appear here (with their original description) for reference.

### [HIGH] (algorithm — OPEN) A source whose line interferes with a strong unmodeled (NORM/background) line must be co-fit with that interferer
- **Where:** `src/FitPeaksForNuclides.cpp` — source/NORM curve construction + clustering (`cluster_gammas_to_rois`), and the caller's source-selection path (`fit_peaks_for_nuclides` entry); also relevant to the `FitNormBkgrndPeaks` NORM source set (`get_norm_sources`).
- **Issue:** Fitting Eu-152 alone on trinitite lets Eu-152's weak 1457 keV line absorb part of K40's strong 1460 keV line; the joint ROI becomes degenerate, the rel-eff curve is skewed, and spurious high-E "peaks" (1347/1363) survive. The general pattern: any requested source with a line within a few FWHM of a strong line NOT in the current model (K40 1460, Tl208 2614, ambient Cs137 661 / Co60 1173·1332, etc.) is at risk. The unit tests are now patched to fit `{K40, Eu152}` jointly, but the *engine* should detect and handle this automatically.
- **Action (design, not yet implemented):** When assembling the fit for a source, scan its gammas against (a) the supplied background spectrum's strong peaks and (b) a known-strong-NORM table (reuse `is_near_strong_norm_gamma`); for any source line within ~N·FWHM of a strong interferer that is NOT already modeled (no background subtraction covering it, and not in the source/NORM curves), auto-add that interferer's nuclide to the fit (or, in supply-background mode, rely on subtraction). Prefer co-fitting (same solve) over a fit-K40-first ordering. Guard so this never splits a real doublet. This is the proper fix that makes order/one-at-a-time/background-mode results converge for interfering-line sources.

### [HIGH] Initial rel-eff eqn_order chosen from matched-peak COUNT ignores number of source activities, causing under-determined solves that silently fall back
- **Where:** `src/FitPeaksForNuclides.cpp:5385-5458`  _(R4:init-rel-eff)_
- **Issue:** manual_input.eqn_order/eqn_form are selected solely from peaks_matched.size() (buckets 1/2/3/4/>4). The actual free-parameter count in RelActCalcManual is num_fit_activities + eqn_order, and the solver enforces (num_fit_activities + num_rel_eff_pars_fit) <= peaks.size() (RelActCalcManual.cpp:1075) and num_peaks >= (eqn_order+1)+(num_isotopes-1) (RelActCalcManual.cpp:5488). With multiple sources each matched peak adds an activity, so e.g. 3 matched peaks split across 2 nuclides with config.initial_manual_relEff_3peak_eqn_order=1 needs (1)+(2)=3 params for 3 peaks -> dof 0; the GA-tuned non-HPGe config sets initial_manual_relEff_2peak_eqn_order=2 (FitPeaksForNuclides.cpp:8198) which for even ONE source needs activities(1)+order(2)=3 > 2 peaks and ALWAYS throws ErrorInitializing.
- **Action:** Make the bucket selection at FitPeaksForNuclides.cpp:5385-5458 aware of the number of distinct matched nuclides and clamp the chosen eqn_order so that (num_distinct_matched_nuclides + eqn_order) <= peaks_matched.size() for non-physical forms (and the analogous physical-model parameter budget), reducing eqn_order or switching form rather than relying on the constructor to throw. Compute num_distinct_matched_nuclides from peak_match_results.used_isotopes (already available) before configuring manual_input. Additionally: (a) when an order must be clamped, emit a fallback_warning instead of silently degrading; (b) re-run the GA for the non-HPGe config with this constraint enforced so it stops tuning parameters that are dead for low-peak/multi-source cases; (c) consider distinguishing ErrorInitializing (parameter-count infeasibility, which should trigger order reduction) from genuine fit failure (which should trigger the physical retry), so the physical-model retry is not silently masking an avoidable under-determination. This also supports the user's physical-earlier and init-rel-eff themes: deciding to go physical sooner should be an explicit, logged choice rather than an accidental side effect of an always-failing non-physical configuration.

### [HIGH] Fallback activity estimate uses inconsistent live-time scaling between matched and unmatched peaks
- **Where:** `src/FitPeaksForNuclides.cpp:5217-5247`  _(R4:init-rel-eff)_
- **Issue:** In estimate_initial_rois_fallback the matched-peak branch computes estimated_activity = peak_area/(best_br*best_eff) with NO live_time division (comment at 5219 says to stay consistent with cluster_gammas_to_rois), but the unmatched branch computes estimated_activity = estimated_peak_area/(best_br*best_eff*live_time) WITH a live_time division (5238). Both activities are appended to the same source_age_and_acts (5247) and passed to one cluster_gammas_to_rois call (5264), which multiplies activity by photon.numPerSecond*rel_eff to get expected COUNTS (cluster_gammas_to_rois 3677).
- **Action:** Remove the spurious *live_time divisor in the unmatched branch at src/FitPeaksForNuclides.cpp:5238 so it reads estimated_activity = estimated_peak_area/(best_br*best_eff), matching the matched-peak branch (5221) and the cluster_gammas_to_rois/RelActManual counts convention. Optionally add a brief test or debug assertion that the two branches produce activities of the same order for the same gamma yield, to prevent regressions.

### [HIGH] Fitness mixes opposing sign conventions: two reward terms summed into a minimized cost
- **Where:** `target/peak_fit_improve/NuclideConfig_GA.cpp; target/peak_fit_improve/PeakFitImprove.cpp; target/peak_fit_improve/FitPeaksForNuclideDev.h:NuclideConfig_GA.cpp:1144-1147; PeakFitImprove.cpp:1774-1776; FitPeaksForNuclideDev.h:64-67`  _(R8a:GA-wrapper)_
- **Issue:** The per-spectrum foreground score is fg = initial_fit_weights.find_weight + final_fit_score.total_weight + candidate_peak_score.score (PeakFitImprove.cpp:1774; identical in NuclideConfig_GA.cpp:1144). openGA MINIMIZES objective1 (calculate_SO_total_fitness returns objective1; do_nuclide_config_ga sets SOGA). But CombinedPeakFitScore's own doc (FitPeaksForNuclideDev.h:64-67) states find_weight is 'higher is better (more correct peaks found)' while total_weight is 'lower is better'. Confirmed in source: InitialFit_GA.cpp:1335 find_weight starts as count of found def-wanted (+1 each), partial credit added, then '-= initial_fit_extra_peak_punishment(0.75)*num_extra' (line 1375). CandidatePeak_GA.cpp:559/565 add +1/+partial per detected expected peak and line 610 subtracts '0.5*num_extra'. So both terms INCREASE when more wanted peaks are found and DECREASE for extra peaks. Adding them into a minimized objective means: finding a correct wanted peak raises the cost (penalized), and emitting an extra/false peak LOWERS the cost (rewarded by -0.75 and -0.5 per extra). Only total_weight (>=0, capped area mismatch) pulls the right way.
- **Action:** Make the foreground objective consistently lower-is-better before it reaches openGA. Concretely, negate the higher-is-better terms: in the per-spectrum fg_score (PeakFitImprove.cpp:1774-1776 and NuclideConfig_GA.cpp:1144-1146) use `-(find_weight + score) + total_weight`, OR negate total_fg in the ga_eval return (PeakFitImprove.cpp:1831) while keeping total_weight and the +weight*total_bg penalty as positive costs. This matches the sibling-GA convention (InitialFit PeakFitImprove.cpp:1643; CandidatePeak CandidatePeak_GA.cpp:1216). Also re-examine the failed-fit/exception penalty (100.0 at PeakFitImprove.cpp:1760,1788): after the sign fix it must remain a strictly worse (larger) cost than any successful fit, so pick a penalty that dominates the best achievable negated score rather than a hardcoded 100. Verify by a short GA run: best_total_cost should trend negative (as the other GAs do, e.g. Best=-19.75) and best genomes should increase found wanted peaks and reduce extras. Fix this BEFORE launching the large per-detector-type runs.

### [HIGH] Failed foreground fit returns fixed cost 100.0 that is NOT comparable to the success-path score scale
- **Where:** `target/peak_fit_improve/PeakFitImprove.cpp:1759-1760, 1786-1789`  _(R8a:GA-wrapper)_
- **Issue:** score_one_spectrum returns {100.0, 0.0} on non-Success status and {100.0,0.0} on exception. But the success-path fg score (find_weight+total_weight+score) is, under the sign bug, often a SMALL or even NEGATIVE number (find_weight/score reward terms reduce a minimized cost; total_weight is bounded). A successful fit that finds several peaks could yield fg well below 100, but a successful fit that finds nothing could also be near 0 - meaning a total fit FAILURE (100) is scored as far worse than a successful-but-useless fit. Per-spectrum 100 penalties also dominate the summed objective and are uncapped relative to the ~O(1-20) per-spectrum success contributions, so a config that fails on a few hard spectra is punished far more than one that quietly mis-fits many.
- **Action:** Fix the objective composition in ga_eval/score_one_spectrum (PeakFitImprove.cpp ~1774-1776, and the mirrored reporting path NuclideConfig_GA.cpp:1144-1146): since the GA minimizes (SOGA), make all components costs - negate/subtract the reward terms find_weight and candidate_peak_score.score (or invert the whole objective). Then recalibrate the failure/exception penalty at PeakFitImprove.cpp:1759-1760 and 1786-1789 to the corrected per-spectrum cost scale instead of the magic 100.0 (e.g. a documented worst-case-per-spectrum constant), and consider normalizing the summed objective by spectrum count so a few hard failures cannot dominate. Re-run the per-resolution GA tuning after the fix, since previously tuned configs were optimized against the inverted objective.

### [MEDIUM] Wrong operand in desperation-trial chi2 comparison: re-tests `solution`, not `trial_solution`
- **Where:** `src/FitPeaksForNuclides.cpp:6816-6818`  _(R1:orchestration)_
- **Issue:** The branch that escalates to a FramPhysicalModel desperation fit is gated by `((trial_solution.m_status != Success) || ((solution.m_chi2/solution.m_dof) > 10.0))`. The second clause uses `solution` (the original energy-cal-on fit) rather than `trial_solution` (the just-computed energy-cal-off fit). So the decision to try Physical Model depends on the chi2 of the WRONG solution. If the original `solution` was bad (>10) but the no-ecal `trial_solution` is now good, the code still fires the Physical-Model desperation attempt needlessly; conversely if `solution` happened to be <10 but `trial_solution` is bad, it will not escalate. The acceptance check at 6870-6872 also compares `trial_solution` against `solution` correctly, so only the escalation gate is wrong. Almost certainly a copy/paste bug from the outer condition at 6797-6799 which legitimately uses `solution`.
- **Action:** Change the second clause of the escalation gate at src/FitPeaksForNuclides.cpp:6817 from `(solution.m_chi2 / solution.m_dof) > 10.0` to `(trial_solution.m_chi2 / trial_solution.m_dof) > 10.0`. While there, consider guarding the division (m_dof can be 0 on a failed/degenerate solve) by gating the chi2 clause behind a Success status check, e.g. `(trial_solution.m_status == Success && trial_solution.m_dof > 0 && (trial_solution.m_chi2/trial_solution.m_dof) > 10.0)`, mirroring the safer `std::max(dof,1)` pattern used in the manual-solution path (lines 5496/5605). Re-run the NuclideConfig_GA per CoarseResolutionType afterward, since this changes when the Physical-Model desperation path triggers and will shift fitness.

### [MEDIUM] Production overload unconditionally enables global debug flag -> data race under GA's parallel workers
- **Where:** `src/FitPeaksForNuclides.cpp:8299-8302, 8593`  _(R1:orchestration)_
- **Issue:** `local_debug_printout` is a plain (non-atomic) file-scope bool (defined line 85). The main fit_peaks_for_nuclides does `const bool prev_debug = local_debug_printout; local_debug_printout = true;` at entry and restores `prev_debug` at 8593, with a TODO saying it is a temporary Xe133 debug hack. eval_peaks_for_nuclide (FitPeaksForNuclideDev.cpp 317-364) and the GA run up to 3 worker threads each calling this function concurrently. Concurrent read-modify-write of `local_debug_printout` is a data race (UB), and the save/restore is itself racy: worker A sets true, worker B saves 'true' as its prev and restores 'true', leaving the flag stuck on after all return. should_debug_print() is consulted in dozens of hot paths.
- **Action:** Remove the temporary debug hack at src/FitPeaksForNuclides.cpp:8301-8302 and the restore at 8593 (revert to the quiet default) before trusting any GA optimization run. If transient tracing must stay, make local_debug_printout thread_local (so each worker has its own copy) AND wrap the set/restore in an RAII guard so every return path (including the early returns at 8310 and 8318) restores the prior value. Do not leave an unsynchronized read-modify-write of a shared non-atomic bool across the openGA worker threads.

### [MEDIUM] FWHM-fit catch block sets got_fwhm_fcn=true and the 'use generic' comment is wrong
- **Where:** `src/FitPeaksForNuclides.cpp:8372-8376`  _(R1:orchestration)_
- **Issue:** On any std::exception from performResolutionFit/removeOutlyingWidthPeaks the catch logs '...WIll use generic FWHM.' but sets `got_fwhm_fcn = true`, which then SKIPS the generic-detector fallback block (8379-8410). If the FIRST performResolutionFit threw (e.g. degenerate widths), fwhm_coefficients is left empty and the function bails out at the empty-check 8422 with FailedToSetupProblem instead of falling back to a generic detector. If the second performResolutionFit threw after the first succeeded, stale first-pass coefficients (pre outlier-removal) are used silently with no warning. Either way the code does NOT do what the comment claims.
- **Action:** Change the catch block at src/FitPeaksForNuclides.cpp:8372-8376 to set got_fwhm_fcn=false so the generic-detector fallback (8379-8410) executes when no usable coefficients were produced. To also handle the second-pass failure cleanly, before deciding, check whether fwhm_coefficients is already non-empty and finite (first pass succeeded): if so, keep those first-pass coefficients, set the energy range from the pre-outlier peaks, and push a local_warnings entry noting outlier removal/refit failed; otherwise fall back to the generic detector. Also fix the "WIll" typo in the log message.

### [MEDIUM] Desperation physical-model escalation gated on solution (not trial_solution) chi2/dof
- **Where:** `src/FitPeaksForNuclides.cpp:6816-6818`  _(R2:engine-A)_
- **Issue:** After running the NoFit-energy-cal retry (trial_solution, line 6810), the inner branch that escalates to the FramPhysicalModel 'desperation' fit checks ((trial_solution.m_status != Success) || ((solution.m_chi2/solution.m_dof) > 10.0)). The chi2 term references the ORIGINAL solution, not trial_solution. Since the outer if (6797-6799) was entered (for a Success-status solution) precisely because solution.m_chi2/solution.m_dof>10.0, this inner condition is essentially always true. So the desperation physical-model solve runs even when the no-ecal trial_solution was excellent. Looks like a copy/paste of the outer predicate.
- **Action:** Change line 6817 from `(solution.m_chi2 / solution.m_dof) > 10.0` to `(trial_solution.m_chi2 / trial_solution.m_dof) > 10.0` so the escalation triggers on the just-computed no-ecal retry quality, matching the comment at 6814 and the intent of the 'physical-earlier' tuning. Optionally, preserve the no-ecal trial_solution before overwriting it at 6867 and have the final selection (6870-6872) choose the best of {original solution, no-ecal trial, desperation result} so a good no-ecal fit cannot be clobbered by a worse desperation fit.

### [MEDIUM] Division by zero / NaN in chi2/dof gating when m_dof == 0
- **Where:** `src/FitPeaksForNuclides.cpp:6799, 6817, 6872`  _(R2:engine-A)_
- **Issue:** RelActCalcAuto sets m_dof to 0 when num_residuals <= num_effective_parameters (RelActCalcAuto.cpp:6698). All three gating expressions compute solution.m_chi2 / solution.m_dof (and trial_solution.m_chi2 / trial_solution.m_dof) with no guard. With dof==0 this is inf or nan. inf>10.0 is true (over-triggers retries), but nan>10.0 is false (skips a retry that should happen) and nan in the final comparison (6872) makes 'trial better than solution' always false, so a genuinely better physical-model/no-ecal solution can be silently rejected.
- **Action:** Add a small helper, e.g. reduced_chi2(sol) returning (sol.m_dof > 0 ? sol.m_chi2 / sol.m_dof : std::numeric_limits<double>::max()) so a dof==0 (degenerate/over-determined) fit is treated as 'bad' deterministically, and use it at FitPeaksForNuclides.cpp:6799, 6817, and 6872. Separately fix line 6817 to evaluate trial_solution (the no-ecal result just computed) instead of the stale original solution. After this, the final selection at 6872 will reliably prefer a successful trial whenever the original is failed or has worse/undefined reduced chi2, instead of silently keeping the worse solution due to nan comparisons. Critical files: /Users/wcjohns/coding/InterSpec_peaks_for_source_opt/src/FitPeaksForNuclides.cpp (lines 6797-6877), /Users/wcjohns/coding/InterSpec_peaks_for_source_opt/src/RelActCalcAuto.cpp (line 6698, m_dof assignment), /Users/wcjohns/coding/InterSpec_peaks_for_source_opt/InterSpec/RelActCalcAuto.h (lines 1456,1462, m_chi2/m_dof types).

### [MEDIUM] ROI-trimming FWHM evaluation uses std::clamp(mean, fwhm_lower, fwhm_upper) without the have_fwhm_range guard used elsewhere
- **Where:** `src/FitPeaksForNuclides.cpp:6246, 6263, 6285-6288, 6451-6453`  _(R2:engine-A)_
- **Issue:** fwhm_lower_energy/fwhm_upper_energy are initialized to -1.0 at the call site (line 8340) and only set from peaks or DRF. They can legitimately remain -1.0 (FWHM fit catch branch at 8372-8376 sets got_fwhm_fcn=true but leaves the energy range untouched; or auto_search_peaks_dq becomes empty after outlier removal at 8365) or be equal (single auto-search peak: 8367-8368 set lower==upper==mean). The four trimming sites call std::clamp(peak_mean, fwhm_lower_energy, fwhm_upper_energy) unconditionally, unlike the merge-ROI block (6707-6712) and other helpers (e.g. 5187-5194) that guard with have_fwhm_range. With lower==upper==-1 the clamp is defined but yields FWHM evaluated at -1 keV (garbage); the margin/abutment math then uses a nonsense FWHM. If the two ever invert (lo>hi) std::clamp is undefined behavior.
- **Action:** Apply the same have_fwhm_range guard used at 6707-6711 to the four trimming sites (6246, 6263, 6285, 6451): compute have_fwhm_range = (fwhm_lower_energy>0 && fwhm_upper_energy>0 && fwhm_lower_energy<fwhm_upper_energy) once, then evaluate peakResolutionFWHM at have_fwhm_range ? std::clamp(energy, lower, upper) : energy. This both removes the bogus -1 keV evaluation and avoids the theoretical inverted-range clamp precondition. Optionally also tighten the catch branch at 8372-8376 so a failed FWHM fit does not silently set got_fwhm_fcn=true while leaving the energy range at -1 (e.g. fall through to the generic-detector branch which sets a valid lower/upper from the DRF).

### [MEDIUM] Desperation physical-model gate re-tests the original solution chi2 instead of the no-ecal trial
- **Where:** `src/FitPeaksForNuclides.cpp:6816-6818`  _(R3:engine-B)_
- **Issue:** After computing trial_solution with energy cal disabled, the nested condition that triggers the FramPhysicalModel desperation solve reads `((trial_solution.m_status != Success) || ((solution.m_chi2 / solution.m_dof) > 10.0))`. The second clause refers to the ORIGINAL `solution`, not `trial_solution`. So if the no-ecal trial succeeded and fits well but the original cal-fit solution was bad (>10), the physical-model desperation solve still fires; conversely a bad trial with a (coincidentally) good original would not. The final acceptance at 6870-6872 then compares trial_solution against solution as expected, so the bug mainly causes an unnecessary extra PhysicalModel solve (cost) or, since trial_solution is overwritten at 6867, can replace a good LnXLnY no-ecal trial with a worse physical-model result that is only kept if it actually beats `solution`.
- **Action:** Change the second clause at line 6817 to reference `trial_solution` rather than `solution`, i.e. gate on `((trial_solution.m_status != Success) || ((trial_solution.m_chi2 / trial_solution.m_dof) > 10.0))`, and guard against `trial_solution.m_dof == 0`/non-success before dividing. For the physical-earlier theme, also promote the hard-coded 10.0 thresholds at 6799 and 6817 into separately tunable NuclideConfigSolution genome parameters so the GA can independently control (a) when to drop energy-cal fitting and (b) when to fall back to the Physical model.

### [MEDIUM] Hard-coded chi2/dof thresholds (10.0) for energy-cal abandon and physical-model fallback are not in the GA genome
- **Where:** `src/FitPeaksForNuclides.cpp:6797-6799, 6816-6817, 7247`  _(R3:engine-B)_
- **Issue:** The decisions to abandon energy-cal fitting and to escalate to a PhysicalModel are gated on the literal `(m_chi2/m_dof) > 10.0` (6799, 6817) and the refinement-accept uses `chi2/dof < solution chi2/dof` (7247). These constants are not members of PeakFitForNuclideConfig (header lists initial_manual_rel_eff_max_chi2_dof and roi_significance_* but no chi2/dof cap for the auto stage). They therefore cannot be optimized per CoarseResolutionType. A value of 10 that is reasonable for HPGe is often unattainable for NaI/CZT where chi2/dof is structurally higher, so NaI fits may always trip the desperation paths (or never, depending on data).
- **Action:** Promote both auto-stage chi2/dof gates to tunable PeakFitForNuclideConfig members (e.g. auto_abandon_ecal_max_chi2_dof and auto_physical_model_escalation_max_chi2_dof) so the GA can optimize them per CoarseResolutionType, mirroring the existing initial_manual_rel_eff_max_chi2_dof. While there, fix the latent bug at :6817 to test trial_solution.m_chi2/trial_solution.m_dof (the no-ecal solve) instead of solution. Leave the relative accept comparison at :7247 as-is (it is not a fixed-threshold issue).

### [MEDIUM] compute_observable_peaks runs on fitted-cal foreground while input peaks are still spectrum-cal (active build)
- **Where:** `src/FitPeaksForNuclides.cpp:7588-7595, 2772-2832`  _(R3:engine-B)_
- **Issue:** OBSERVABLE_PEAKS_USING_ORIGINAL_CAL_WITH_BACK_SUB is 0 (line 30), so the active path is 7589-7595: compute_observable_peaks(result.fit_peaks, foreground, ...) is called with the energy-cal-ADJUSTED foreground (reassigned at 6947), but result.fit_peaks are still in spectrum/original cal (built from m_peaks_without_back_sub) and are only translated to orig cal AFTER (7592-7594). Inside compute_observable_peaks the initial significance filter integrates data at foreground->gamma_integral(mean-fwhm, mean+fwhm) with mean in the wrong cal (2809), and refitPeaksThatShareROI_LM refits spectrum-cal peaks against fitted-cal data (2952). When apply_energy_cal_between fit a non-trivial shift, the data window is offset from the true peak, biasing significance and the refit.
- **Action:** Make the data and peaks use the same cal in the active observable-peaks path. Simplest fix that matches the existing intent: in the #else branch translate result.fit_peaks (and uncombined_fit_peaks) to original cal BEFORE calling compute_observable_peaks, and pass orig_foreground instead of foreground (i.e., behave like the #if branch but without background subtraction), then translate only the resulting observable_peaks if any further adjustment is needed. Alternatively, before line 7590 translate the input peaks from orig cal to fitted cal so they match the fitted-cal foreground, run compute_observable_peaks, then translate fit_peaks/uncombined/observable back. Either way the gamma_integral significance window (2809) and refitPeaksThatShareROI_LM data (2952) must see peak means in the same cal as the Measurement. Add a PERFORM_DEVELOPER_CHECKS assert that, when apply_energy_cal_between is active, the peak means handed to compute_observable_peaks are consistent with foreground->energy_calibration(). Confirm by running the NuclideConfigGA HPGe/NaI cases with fit_energy_cal=true vs DoNotVaryEnergyCal and checking that observable-peak counts/areas no longer diverge purely due to the cal toggle.

### [MEDIUM] Final residual energy-cal adjustment from the last refine solve is never propagated to translation
- **Where:** `src/FitPeaksForNuclides.cpp:6905-6965, 7310-7311, 7554-7576`  _(R3:engine-B)_
- **Issue:** Each refine iteration applies the PRIOR solution's fitted cal to `foreground` (6907,6947), then re-solves on that already-adjusted foreground (7311). The new solution.m_energy_cal_adjustments are residuals relative to that adjusted foreground, but the loop only re-applies cal at the TOP of the next iteration. On the final accepted iteration there is no subsequent top-of-loop, so the last solve's residual shift is folded into m_peaks_without_back_sub (spectrum cal = last input foreground cal) yet the translate step (7560-7561) maps from foreground->energy_calibration() (the cal set at the START of that iteration) to orig_cal. The final residual is thus dropped from the translation, leaving peak means slightly off in the returned original-cal peaks.
- **Action:** In translate_peaks_to_orig_cal (FitPeaksForNuclides.cpp:7560) source the `fitted_cal` from solution.m_foreground->energy_calibration() rather than the loop-local `foreground` variable, so the translation source cal always matches the cal the returned peaks (solution.m_peaks_without_back_sub) actually live in. Alternatively, after the refine loop, reset `foreground = solution.m_foreground` (or only advance `foreground` at the END of an iteration once the solution is accepted, rather than at the top using the prior solution). A PERFORM_DEVELOPER_CHECKS assert that *foreground->energy_calibration() == *solution.m_foreground->energy_calibration() at the translate site would catch regressions.

### [MEDIUM] combine_overlapping_peaks_in_rois merges distinct-source peaks and keeps only the dominant peak's source
- **Where:** `src/FitPeaksForNuclides.cpp:2434-2490, 2546-2564`  _(R3:engine-B)_
- **Issue:** Greedy clustering combines a smaller peak into a larger anchor whenever should_combine_peaks returns true (within 1.5*sigma OR smaller_area < 4*sqrt(tail contribution), 2407-2428). combine_peaks copies the dominant (largest-area) peak (2472), so the merged peak carries ONLY the dominant source attribution; the smaller peak's source line is lost. should_combine_peaks does not consider whether the two peaks belong to different sources. Near clustered multiplets (and around the ~120 keV knee where many lines pile up) a genuine distinct source line riding the tail of a larger neighbor can be absorbed and mis-attributed.
- **Action:** Make merging source-aware. Two concrete options: (a) In should_combine_peaks / the clustering loop (src/FitPeaksForNuclides.cpp:2546-2564), refuse to combine when the candidate peak has a different source assignment than the anchor AND the candidate is itself statistically significant (e.g. its own area exceeds a detection threshold), so only true tail-contamination peaks (or same-source lines) get absorbed. (b) If a merge of differing-source peaks does occur, do not silently keep only the dominant source: either record the secondary source(s) on the merged PeakDef (the header at InterSpec/RelActCalcAuto.h:110-112 already calls for a multi-source-per-peak representation) or keep the smaller distinct-source peak separate for attribution/scoring purposes. At minimum, only allow cross-source merging under criterion 2 (overwhelmed-by-tail), never under the unconditional 1.5*sigma criterion 1 (line 2411), since criterion 1 will merge a fully significant neighboring line of another source purely on proximity. Given the GA scores observable_peaks per source, this directly affects the knee/multiplet regions (~120 keV) the user is concerned about.

### [MEDIUM] rel-eff clamping flattens the curve outside the matched-peak energy span, mispredicting ROI counts for far gammas
- **Where:** `src/FitPeaksForNuclides.cpp:5627-5635`  _(R4:init-rel-eff)_
- **Issue:** manual_rel_eff clamps energy to [peaks_matched.front().m_energy, peaks_matched.back().m_energy] and returns the endpoint rel-eff outside that range. cluster_gammas_to_rois (3673) requests rel_eff at ALL source gamma energies, including gammas well below/above the matched-peak span. Clamping yields a constant (often the low-energy-high value) for high-energy gammas where the true efficiency is much lower.
- **Action:** Replace flat endpoint clamping with a physically-motivated extrapolation outside the matched span instead of a constant. Options, in rough order of robustness: (1) When extrapolating, fall back to the DRF intrinsic-efficiency shape (or the FramPhysical/Empirical model curve) anchored to the endpoint rel-eff value, so high-energy gammas get a declining efficiency rather than the flat in-span minimum — this also ties into the team's 'physical-earlier' theme. (2) At minimum, gate ROI creation in cluster_gammas_to_rois on whether the gamma energy is within (or close to) the matched-peak span: for clamped/extrapolated gammas, require a stronger passes_data_area / passes_signif margin before seeding a ROI, so a purely-predicted far-energy gamma cannot create a ROI on prediction alone. (3) Track per-gamma an 'extrapolated' flag and exclude extrapolated far gammas from the spurious-ROI-prone acceptance entirely, letting them only ride along inside ROIs anchored by in-span gammas. Add a targeted regression test for a narrow-span single-source case (e.g. one matched peak ~120 keV with the source also emitting a high-energy gamma) to confirm no spurious high-energy ROI is seeded.

### [MEDIUM] chi2/dof retry trigger and acceptance are unreliable when dof <= 0
- **Where:** `src/FitPeaksForNuclides.cpp:5496-5614`  _(R4:init-rel-eff)_
- **Issue:** chi2_dof = m_chi2/max(m_dof,1) and final_chi2_dof = m_chi2/max(m_dof,1). m_dof from RelActCalcManual.cpp:5518 is num_peaks-(eqn_order+1)-(num_isotopes-1) and can be 0 or negative for the order/source combos this code produces; max(dof,1) hides that. So a fit with zero true degrees of freedom (interpolating) can show an artificially small chi2/dof and pass both the >20 retry gate and the initial_manual_rel_eff_max_chi2_dof (default 25, non-HPGe tuned) acceptance, OR a marginally-determined fit can be compared against the physical retry on a non-comparable chi2 scale (different residual normalization, USE_RESIDUAL_TO_BREAK_DEGENERACY adds an extra residual for non-physical only, number_residuals 1080-1089).
- **Action:** Treat m_dof <= 0 explicitly instead of hiding it with max(m_dof,1). Options: (1) when manual_solution.m_dof <= 0, do not trust chi2/dof as a goodness gate -- either force the physical-model retry, or reject the fit as degenerate (require strictly more matched peaks than fit parameters, e.g. require m_dof >= 1 for the source-presence acceptance), so an interpolating fit cannot be accepted as if absent sources were present. (2) For the retry-vs-original comparison at line 5584, compare on raw m_chi2 (or a properly degenerate-safe reduced chi2) since the numerators are comparable, and only fall back to dof-normalized comparison when both dof are positive. (3) Surface m_dof in the debug output (already printed) into the decision so the GA tuning of initial_manual_rel_eff_max_chi2_dof is not optimizing against masked zero-dof points.

### [MEDIUM] data_area / significance use gross counts, not continuum-subtracted, so a strong nearby (other-source or NORM) peak inflates significance and rescues an otherwise-empty ROI
- **Where:** `src/FitPeaksForNuclides.cpp:3779-3797`  _(R5:clustering)_
- **Issue:** data_area = foreground->gamma_integral(lower,upper) is the gross integral over the window, and signif = counts_in_region/sqrt(data_area). counts_in_region only sums the source's expected counts, but data_area is dominated by whatever real counts (continuum + unrelated peaks) sit in the window. A window that happens to overlap a large unrelated peak (e.g. a NORM 1460 keV, or another source's strong line) gets a LARGE data_area, which (a) easily passes passes_data_area and (b) for a fixed counts_in_region, REDUCES signif (counts/sqrt(area)), so it is not a false-positive rescue but the opposite: a genuine weak source line sitting on a Compton shoulder or near a big neighbor gets its significance suppressed and the ROI dropped. The min_data_area_keep gate also passes trivially anywhere the spectrum has counts, so it provides almost no discrimination. Net: the significance metric conflates 'is there signal here' with 'is the background high here', biasing rejection in busy/high-background regions exactly where background-present vs background-absent behavior should differ.
- **Action:** Replace the gross data_area in the significance denominator with a continuum-subtracted estimate. Concretely, in cluster_gammas_to_rois (around 3779-3795) estimate the local continuum under [lower,upper] using the existing left/right sideband approach (mirroring should_use_step_continuum at 3494-3521), subtract it to get an estimated net/excess area, and compute signif against sqrt(continuum-or-total) Poisson noise so the metric measures signal-over-noise rather than signal-over-gross-counts. Apply the same fix to the duplicate metric at 4511-4513. This removes the background-dependence of the keep/drop decision (directly serving the background-invariant-fitting goal) and stops genuine weak lines near strong NORM/other-source peaks from being suppressed. Keep the min_data_area_keep gate only as a near-empty-window guard. Re-tune the three keep thresholds via the NuclideConfig GA after the denominator semantics change, since their current GA-tuned values (e.g. 8208-8223) were optimized against the old gross-area definition.

### [MEDIUM] Single shared cluster_num_sigma / roi_width / FWHM treatment across the whole energy range — no knee, so low-energy (steep rel-eff, x-ray dense) and high-energy regions are clustered with identical rules
- **Where:** `src/FitPeaksForNuclides.cpp:3756-3757, 3907-3908, 4172`  _(R5:clustering)_
- **Issue:** Window half-width (cluster_num_sigma*sigma, 3756-3757), ROI extent (roi_width_num_fwhm_*, 3907-3908) and break width (max_fwhm_width*mid_fwhm, 4172) all scale only with the local FWHM and fixed global multipliers. Below ~120 keV, lines are dense (x-rays, escape, backscatter, U/Pu complexes), the rel-eff curve is steep and least-trustworthy, and FWHM is small so windows are narrow in keV but lines are tightly packed; above ~120 keV lines are sparser and rel-eff smoother. The same num_sigma therefore over-merges below the knee and the same significance thresholds are applied where the rel-eff estimate (used to compute counts_in_region at 3673,3677) is least reliable. There is no mechanism to (a) split the energy axis at a chosen knee and cluster each segment with its own settings, nor (b) let the rel_eff_fcn be a piecewise/knee curve.
- **Action:** Treat as the anchor finding for the 'knee' theme rather than a bug fix. Two complementary options, both implementable here: (1) Make the clustering policy energy-aware — add a knee energy (GA-tunable per CoarseResolutionType) and let GammaClusteringSettings carry below-knee vs above-knee variants of cluster_num_sigma / roi_width_num_fwhm_* / significance thresholds, selecting per-gamma in cluster_gammas_to_rois at 3756-3757, 3907-3908, 4172, 3793-3797. This directly counters over-merging and unreliable rel-eff-gated acceptance below ~120 keV. (2) Allow rel_eff_fcn to be piecewise/knee-shaped (a discontinuity or separate segment fit at the knee) so the amplitudes feeding 3673-3677 -> counts_in_region are more reliable at low energy. Option (1) is the lower-risk, higher-leverage first step because it is local to this function and naturally exposed to the existing NuclideConfigSolution genome. Also worth coupling with the 'init-rel-eff' and 'fwhm-constraint' themes since both feed the same acceptance gate.

### [MEDIUM] should_use_step_continuum keys off a single max-amplitude gamma and fixed +/-1.5 FWHM windows that can land on neighbor peaks
- **Where:** `src/FitPeaksForNuclides.cpp:3475-3524`  _(R6:observable-peaks)_
- **Issue:** It picks the single largest-amplitude gamma (3476-3478) and integrates raw foreground in 0.25-FWHM windows centered at ref +/- 1.5 FWHM (3496-3513). In a multi-gamma cluster (which is exactly when ROIs are wide), the +/-1.5 FWHM window can fall on an adjacent real peak's tail or on another gamma's photopeak, so left_sum/right_sum reflect a neighbor rather than the continuum step. The nsigma test (3521) then mis-fires. Also if either window extends beyond the ROI it bails returning false (3506-3507), so legitimately stepped narrow ROIs (peak near ROI edge) never get a step continuum.
- **Action:** Redesign the step-continuum decision so it is robust and background-aware: (a) Decouple the probe geometry from the ROI/cluster width or scale the cluster width when a step test is wanted; at minimum, when a probe window would exceed the ROI, clamp/shrink the window or measure the step just outside the ROI edges instead of returning false at 3506-3507, so strong isolated peaks (e.g. Cs-137 662) can be assigned a step. (b) Place the left/right probe windows relative to the ROI edges or relative to the gap between the dominant gamma and its nearest in-ROI neighbor, and explicitly skip/avoid windows that overlap any other cluster gamma's +/- ~1.5 sigma region (the other gamma positions are already available in cluster.gamma_energies) to prevent neighbor-peak contamination. (c) Make the decision background-consistent: subtract the provided background spectrum (when present) before integrating, or base the step test on the fitted/estimated continuum rather than raw foreground counts, so the with/without-background results agree. (d) Consider letting the subsequent RelActCalcAuto fit choose between smooth and step continuum families (or compare chi-square of both) instead of a fixed +/-1.5 FWHM heuristic, since the continuum type materially affects fitted peak area. As this is a read-only review, no code was changed; these are recommendations for the implementation/GA-tuning phase. Also worth exposing the 1.5/0.25 FWHM probe offsets as tunable NuclideConfigSolution genome parameters so the GA can correct the geometry per CoarseResolutionType.

### [MEDIUM] Escape-peak ROIs sized with parent FWHM, not escape-energy FWHM, and added only in the no-background path
- **Where:** `src/FitPeaksForNuclides.cpp:1274-1373`  _(R7:roi-norm-escape)_
- **Issue:** roi_width_lower/upper = config.auto_rel_eff_roi_width_num_fwhm_* * parent_peak->fwhm(), and the per-side minimum is 0.5*parent_peak->fwhm(). For Th232 2614 keV the SE (2103) and especially DE (1592) lines are ~1000 keV lower where HPGe FWHM is substantially smaller, so the SE/DE ROIs are built too wide. The added RoiRange leaves continuum_type=Quadratic and range_limits_type=CanBeBrokenUp at defaults. Because add_escape_peak_floating_peaks_if_appropriate() returns early unless fit_norm_peaks (line 1136), these ROIs/floating peaks exist ONLY when NORM peaks are fit, i.e. only in the background-absent (NoBackgroundFitNorm) GA mode.
- **Action:** Size the SE/DE escape ROIs (and their 0.5-FWHM minimum) using the detector FWHM evaluated at the escape energy (se_energy/de_energy) rather than parent_peak->fwhm() - e.g. derive expected FWHM from the DRF or the FWHM functional fit at those energies, mirroring how the main ROI path uses each edge peak's own FWHM at FitPeaksForNuclides.cpp:2737-2740. Separately, reconsider the fit_norm_peaks gate at FitPeaksForNuclides.cpp:1136: escape peaks of high-energy source lines (Th232 2614) are a property of the source, not of NORM fitting, so adding the escape ROIs/floating peaks should not depend on whether NORM peaks are being fit; decouple it so background-present and background-absent runs treat the 1592/2103 keV region consistently. Lower priority: set an explicit continuum_type (e.g. Linear for these narrow escape ROIs) instead of relying on the Quadratic default.

### [MEDIUM] resolve_overlapping_rois only extends upper_energy; drops the merged ROI's continuum_type/range fields and never re-splits
- **Where:** `src/FitPeaksForNuclides.cpp:1467-1534`  _(R7:roi-norm-escape)_
- **Issue:** On overlap it does last.upper_energy=max(last.upper_energy,current.upper_energy) and discards `current` entirely. The merged ROI keeps `last`'s continuum_type (e.g. Linear/step from a source ROI) and silently drops `current`'s continuum_type — if the escape/source ROIs had different continuum settings, the merge result depends on sort order (ascending lower_energy) rather than on which ROI is dominant. There is no width cap on the merge (contradicting the doc comment 'merging them if the combined width is reasonable'), so a chain of overlaps can produce one very wide ROI. The PERFORM_DEVELOPER_CHECKS branch only asserts overlaps are escape-related, masking other ROI-generation overlaps in release builds.
- **Action:** When merging, explicitly decide which ROI's continuum_type and range_limits_type the result should take (prefer the source/Fixed ROI over a default-valued escape ROI rather than relying on lower_energy sort order), e.g. pick the more conservative/dominant continuum and keep range_limits_type=Fixed if either input was Fixed. Set continuum_type/range_limits_type on the escape ROIs at 1290-1366 instead of leaving them at struct defaults. Add the width cap the doc promises (and either implement the split path or fix the doc comment to say it only merges). Optionally promote the non-escape-overlap warning out of PERFORM_DEVELOPER_CHECKS so unexpected overlaps are visible in release builds.

### [MEDIUM] openGA RNG is seeded from wall-clock; runs are non-deterministic and unseeded by CLI
- **Where:** `target/peak_fit_improve/build_vscode/_deps/openga-src/src/openGA.hpp; target/peak_fit_improve/NuclideConfig_GA.cpp:openGA.hpp:371-374; NuclideConfig_GA.cpp:1556-1576`  _(R8a:GA-wrapper)_
- **Issue:** openGA's ctor seeds rng from std::chrono::high_resolution_clock (openGA.hpp:372-374). do_nuclide_config_ga never sets a seed and there is no CLI option to fix one (PeakFitImprove.cpp GA options around 1026-1058 expose population/generations/threshold but no seed). The mutate() fallback after 1000 tries also uses std::random_device (NuclideConfig_GA.cpp:753).
- **Action:** Add a `--ga-seed` CLI option (boost::program_options) in PeakFitImprove.cpp and thread it into do_nuclide_config_ga so the GA's rng is reseeded deterministically right after `GA_Type ga_obj;` is constructed and before solve() (NuclideConfig_GA.cpp ~1556). Since openGA exposes no setter, either set ga_obj.rng.seed(seed) directly (the member is accessible) or add a thin reseed call. Replace the std::random_device fallback at NuclideConfig_GA.cpp:753 with a draw from the GA's own rnd01/rng stream so the fallback path stays reproducible. Log the chosen seed to nuclide_config_ga_results.txt / stdout so any run can be reproduced. For full determinism also document that bit-identical replay additionally requires single-threaded eval (N_threads=1) or a per-individual RNG-draw scheme, since multi_threading=true can reorder rng consumption. Default to a random seed (preserving current behavior) but always print it.

### [LOW] Misleading 'best RelEff curve type selected' messaging implies a ranking that does not exist
- **Where:** `target/peak_fit_improve/FitPeaksForNuclideDev.cpp:417-418, 482-483`  _(R1:orchestration)_
- **Issue:** The dev driver prints 'Best RelEff curve type solution (...) selected with chi2/dof=' and the failure path says 'fit_peaks_for_nuclides failed for all RelEff curve types', but fit_peaks_for_nuclides only ever fits the single config.rel_eff_eqn_type (5770) plus opportunistic Physical-Model desperation retries. No set of curve types is enumerated or chi2-ranked. The header doc (FitPeaksForNuclides.h, function comment) similarly says step 5 'Calls fit_peaks_for_nuclide_relactauto with three different configurations / Returns the best result based on chi2', which does not match the implementation.
- **Action:** Treat as a documentation/observability fix, not a correctness bug. (1) Update the two dev-driver strings (FitPeaksForNuclideDev.cpp:417-418 and 482-483) to say the single configured curve type was used, e.g. 'RelEff fit (<name>) chi2/dof=...' and 'fit_peaks_for_nuclides failed (curve type <name>)', or keep wording only if a real multi-curve comparison is implemented. (2) Fix the header doc FitPeaksForNuclides.h:407-408 to describe the actual single-config + conditional Physical-Model desperation flow. (3) Separately, for the 'physical-earlier' theme the team should decide whether to add a genuine curve-type comparison (try LnXLnY/Physical and pick best by chi2/dof) rather than relying solely on the GA genome + narrow desperation retries — but that is a design enhancement, tracked apart from this messaging fix.

### [LOW] FitNormBkgrndPeaks forces background subtraction off, coupling NORM-fit and background-present behavior
- **Where:** `src/FitPeaksForNuclides.cpp:8494-8519`  _(R1:orchestration)_
- **Issue:** When FitNormBkgrndPeaks or FitNormBkgrndPeaksDontUse is set, `long_background = nullptr` so background subtraction is disabled for the whole fit, and NORM ROIs are built on a second rel-eff curve. The choice of whether to fit NORM peaks is thus entangled with whether background subtraction happens. The GA background-penalty pass (NuclideConfig_GA.cpp 300-311) sets FitNormBkgrndPeaks for NoBackgroundFitNorm and always passes background=nullptr, so for that pass it is consistent, but in the GUI/foreground path enabling NORM fitting silently turns off background subtraction even if a background spectrum was supplied.
- **Action:** Keep the behavior for the active-NORM-fit case but make the coupling explicit and reconsider the DontUse case. Concretely: (1) Document on the FitNormBkgrndPeaks/FitNormBkgrndPeaksDontUse enum in InterSpec/FitPeaksForNuclides.h:392-397 that setting either flag disables background subtraction for the whole fit. (2) Add a result.warnings entry at src/FitPeaksForNuclides.cpp:8497 when a non-null long_background is being discarded, so the caller/user knows subtraction was turned off. (3) Re-evaluate whether FitNormBkgrndPeaksDontUse should null long_background at all: since those NORM peaks are not returned, the source ROIs could still benefit from background subtraction; consider only nulling for the true FitNormBkgrndPeaks (return-NORM) case. This is orthogonal to the four optimization themes (knee/physical-earlier/init-rel-eff/fwhm-constraint) and should be treated as a low-priority clarity fix, not a fit-quality change.

### [LOW] fwhm_estimation_method hard-coded to StartFromDetEffOrPeaksInSpectrum; not tunable and ignores config
- **Where:** `src/FitPeaksForNuclides.cpp:5948`  _(R2:engine-A)_
- **Issue:** options.fwhm_estimation_method is always StartFromDetEffOrPeaksInSpectrum regardless of detector type or whether a trustworthy DRF FWHM exists. The GA tunes config.fwhm_form but cannot influence the estimation method, and there is no path to Fixed*/StartingFromDetectorEfficiency. For HPGe with a good DRF, fixing FWHM would both stabilize the fit and avoid the unbounded-coefficient concern; for NaI with few peaks, starting purely from peaks can be unstable.
- **Action:** Add an fwhm_estimation_method field to the NuclideConfig struct (InterSpec/FitPeaksForNuclides.h, alongside fwhm_form at line 193) and use it at src/FitPeaksForNuclides.cpp:5948 instead of the hard-coded literal. Add a corresponding gene to the GA genome (NuclideConfig_GA.cpp serialize/init/mutate) constrained per detector type, e.g. allow Fixed*/StartingFromDetectorEfficiency for HPGe with a trustworthy DRF (stabilizes FWHM and addresses the unbounded-coefficient concern) and keep StartFromDetEffOrPeaksInSpectrum for sparse-peak NaI. Note FwhmForm::NotApplicable requires FixedToDetectorEfficiency (per RelActCalcAuto.h:425-430), so the GA must keep these two coupled to avoid invalid combinations. Also guard against StartingFromDetectorEfficiency/FixedToDetectorEfficiency when the DRF lacks resolution info, since those throw.

### [LOW] rel_eff_eqn_order not bounded against number of fit-able lines at the auto stage (only the manual seed is reduced)
- **Where:** `src/FitPeaksForNuclides.cpp:5770-5777`  _(R2:engine-A)_
- **Issue:** config.rel_eff_eqn_order (GA range 0..6) is copied verbatim into rel_eff_curve.rel_eff_eqn_order for non-physical forms with no check against how many independent source lines actually land in ROIs. The solver reduces only the MANUAL seed order when underdetermined (RelActCalcAuto.cpp:3884-3890) and then re-expands to the full order padding higher coefficients with 0 (4078-4088); the non-linear auto fit still carries the full order. A high-order LnXLnY with 1-2 lines is underdetermined and the rel-eff curve can swing wildly between/at line energies, which is the classic cause of mis-attributed amplitudes. Order 0 for a multi-line source is the opposite failure (flat curve).
- **Action:** Add a bound on rel_eff_eqn_order at the auto-fit setup in FitPeaksForNuclides.cpp (around 5770-5777), mirroring the manual-seed logic at RelActCalcAuto.cpp:3884-3890: count the independent fit-able lines/peaks expected to land in ROIs across synced_sources and clamp rel_eff_curve.rel_eff_eqn_order so num_free_pars stays >= some small margin (e.g. order <= max(0, num_independent_lines - num_isos - 1), or simply order < number of distinct line energies). This both prevents an over-flexible curve for sparse sources and reduces GA fitness noise from the same genome order behaving very differently across spectra. Consider also providing a sensible floor (>=1) for clearly multi-line sources, and re-enabling a guard analogous to the commented-out check at RelActCalcAuto.cpp:2940-2941. This dovetails with the user's init-rel-eff and physical-earlier themes: when lines are too few to support any polynomial order, fall back to the Physical model.

### [LOW] Insignificant-ROI peak filter compares spectrum-cal peak continuum bounds against true-energy ROI ranges
- **Where:** `src/FitPeaksForNuclides.cpp:7414-7493`  _(R3:engine-B)_
- **Issue:** insignificant_roi_ranges is built from solution.m_final_roi_ranges (RelActCalcAuto.h:1344 documents these as TRUE energy / fitted cal). The keep/drop loop iterates solution.m_peaks_without_back_sub, whose continuum bounds are in the SPECTRUM cal (RelActCalcAuto.cpp:6181 derives them from m_fit_peaks_in_spectrums_cal). The match at 7487-7488 uses a fixed 1.0 keV tolerance on peak_roi_lower/upper vs roi_range.first/second. When energy cal is fit (sm_energy_cal_multiple=100, RelActCalcAuto.h:1416), the two cals can differ by more than 1 keV (common for NaI/CZT), so peaks that ARE in an insignificant ROI fail the bound match and are NOT flagged as in_insignificant_roi; they then pass to the else-branch where mean_in_roi is checked against m_final_roi_ranges_in_spectrum_cal (correctly spectrum cal) and get KEPT.
- **Action:** Fix is cheap and worthwhile for robustness even though impact is low. Build insignificant_roi_ranges from solution.m_final_roi_ranges_in_spectrum_cal instead of m_final_roi_ranges (FitPeaksForNuclides.cpp:7414-7418) so the 7487-7488 bound match compares spectrum-cal to spectrum-cal; this also removes the dead-code smell. Separately (and more meaningfully for correctness near thresholds), make compute_roi_chi2_significance cal-consistent: pass the spectrum-cal ROI range (m_final_roi_ranges_in_spectrum_cal) when selecting peaks (line 2194) and when computing the channel window (lines 2202-2203), or pass the spectrum used for fitting in the matching cal. Do not raise this as a high-severity returned-peaks bug; verify with a NaI/CZT case where energy cal is fit and a >1 keV shift occurs, confirming the returned peak set is unchanged (which supports the low-severity assessment) while the significance flagging may shift.

### [LOW] Edge-ROI restoration uses front ROI width and BR window heuristics that can misbehave for one-sided spectra
- **Where:** `src/FitPeaksForNuclides.cpp:7044-7179`  _(R3:engine-B)_
- **Issue:** The dropped-edge-ROI restore computes nearest_max_br over a window `half_width = 0.25*(refined_rois.front().upper-front().lower+1.0)` (7134) -- i.e. derived from the FIRST refined ROI's width even when evaluating a high-energy edge ROI restoration whose nearest kept ROI is elsewhere; the window for finding the comparison BR is unrelated to the nearest ROI's own width. The area/BR sanity gate (area_ratio > br_ratio*10, 7154) silently skips when nearest_max_br==0 (7149), allowing restoration with no NORM guard. best_peak_area threshold of 100 (7082) is an absolute count with no live-time/continuum scaling.
- **Action:** Make the edge-ROI restore guard robust rather than front-ROI-anchored. Specifically: (a) size the BR-comparison window from the nearest kept ROI's own width/center (not refined_rois.front()) at line 7134; (b) decide what happens when nearest_max_br==0 at line 7149 -- treat an unfound comparison BR as a reason to be MORE conservative (skip/keep-dropped) rather than silently bypassing the NORM guard, or widen the search to the source's strongest in-window line; (c) scale the best_peak_area threshold at line 7082 by live-time and/or local continuum (or expose it as a CoarseResolutionType-specific GA-tunable in NuclideConfigSolution rather than a hardcoded 100.0). All three are low-priority since this only affects extreme-energy peak survival; (a)+(b) together are the meaningful fix. No correctness invariant is violated, so this can be batched with other rel-eff-extrapolation hardening rather than fixed urgently.

### [LOW] ROI significance veto uses Quadratic continuum-only null even for ROIs fit with step/linear continua
- **Where:** `src/FitPeaksForNuclides.cpp:2252-2278, 2319-2340`  _(R3:engine-B)_
- **Issue:** chi2_with_peaks is fit with the ROI's actual continuum type (cont_type, 2235-2250) but chi2_continuum_only forces Quadratic (2257). For ROIs whose real continuum is a step or higher-flexibility form, the quadratic null is less flexible than the with-peaks continuum, inflating chi2_reduction and making such ROIs look more significant; conversely the separate quad_cont veto (2319-2336) uses max(cont_type, Quadratic). The chi2_with_peaks holds peak amplitudes FIXED (from the fitted, possibly bg-subtracted solution) while re-fitting only the continuum on the raw (cal-adjusted) foreground; a peak amplitude fit on bg-subtracted data but evaluated on non-bg-subtracted data can mis-state the reduction near strong background lines.
- **Action:** Make the chi2-reduction null use the same continuum-type policy as the quad-cont veto, i.e. max(cont_type, Quadratic), so the with-peaks and null fits differ only by the presence of peaks (not by continuum flexibility) and the two significance tests stay consistent; for step-continuum ROIs prefer comparing against the same step type so the step basis function does not by itself inflate the reduction. To remove the with/without-background discrepancy, run the significance comparison on a single consistent dataset/amplitude basis (either evaluate against the background-subtracted spectrum the amplitudes were fit on, or re-fit peak amplitudes on the raw foreground for this null-vs-peaks comparison) rather than mixing bg-subtracted amplitudes with a raw-foreground continuum re-fit. Given low severity, scope this as a hardening change validated through the NuclideConfigGA harness so the GA-tuned thresholds are re-fit after the change.

### [LOW] base_rel_eff_uncert + zero amplitude uncertainty can hard-fail the manual solve
- **Where:** `src/FitPeaksForNuclides.cpp:5312-5320`  _(R4:init-rel-eff)_
- **Issue:** peak_info.m_counts_uncert = peak->amplitudeUncert() and m_base_rel_eff_uncert = config.rel_eff_manual_base_rel_eff_uncert (0.0 for the HPGe header default, 0.308957 for non-HPGe). RelActCalcManual validates (RelActCalcManual.cpp:804-806): if m_counts_uncert==0.0 and m_base_rel_eff_uncert != -1.0 it throws; it also requires base_rel_eff_uncert to be consistently a value in [0,1] or consistently -1 across all peaks (810-816). auto_search_peaks whose amplitudeUncert() is 0 (or not computed) will trip this.
- **Action:** In the conversion loop at src/FitPeaksForNuclides.cpp:5306-5320, guard the per-peak uncertainty before pushing into rel_act_manual_peaks. Concretely: if peak->amplitudeUncert() <= 0.0, either (a) skip the peak from the rel-eff seed, or (b) substitute a sane positive value (e.g. sqrt(amplitude) counting-statistics estimate) so m_counts_uncert stays > 0.0, OR (c) set m_base_rel_eff_uncert = -1.0 consistently across ALL peaks to do an unweighted fit when any peak lacks a usable uncertainty. Whatever is chosen must be applied uniformly so the 808-816 consistency check is not violated. This converts a hard fallback into a usable rel-eff seed. Also reconsider the HPGe default of 0.0 (pure counting-statistics weighting) only as a separate, lower-priority tuning question, not as a correctness fix.

### [LOW] estimate_initial_rois_without_peaks calls merge_rois without foreground or unfit peaks, weakening split/merge quality
- **Where:** `src/FitPeaksForNuclides.cpp:5097`  _(R4:init-rel-eff)_
- **Issue:** merge_rois(initial_rois, config) is called with default empty unfit_auto_peaks and null foreground. Consequently the unfit-peak-between guard (4720-4768) is skipped, and find_spectrum_valley (4808) falls back to the plain midpoint because foreground is null (find_spectrum_valley 3557). The tail-contribution merge check also runs with amplitudes = br*eff scores (gamma.br_times_eff stored as estimated_amplitude at 5087), which are on a different scale than real peak areas used elsewhere.
- **Action:** Add foreground (and the unfit/auto-search peaks the caller already holds) as parameters to estimate_initial_rois_without_peaks and forward them into merge_rois at line 5097, matching the other two merge_rois call sites (6740, 8571). This restores valley-based split-point selection and the interfering-unfit-peak guard in the no-matched-peaks branch. Treat the br*eff-vs-counts amplitude-scale concern as minor since the tail check is ratio-based; no scale fix is required for correctness.

### [LOW] Cluster bounds re-derivation can shrink a window below cluster_num_sigma and silently drop the very gamma that seeded it
- **Where:** `src/FitPeaksForNuclides.cpp:3845-3918`  _(R5:clustering)_
- **Issue:** The keep/reject decision (significance, counts, data_area) is computed on the wide initial window [energy +/- cluster_num_sigma*sigma] (3756-3797). After acceptance, the bounds are REPLACED by weighted_mean +/- roi_width_num_fwhm_{lower,upper}*effective_fwhm (3907-3908). For a cluster of one dominant gamma plus weak satellites, the amplitude-weighted mean is pulled toward the dominant line, but effective_fwhm is inflated by sum_weighted_sq_dev (the spread of means). If roi_width_num_fwhm_lower/upper (manual default 3.0, auto 3.5) is smaller than cluster_num_sigma plus the mean spread, the new ROI can be narrower OR shifted so that a gamma counted in counts_in_region now lies outside [cluster.lower,cluster.upper]. The significance that justified keeping the ROI was computed over a different (often wider) energy span than the ROI actually emitted, so a ROI can be kept on the strength of counts that the final ROI does not even cover.
- **Action:** After re-deriving bounds at 3907-3908, clamp them to still contain the cluster's gamma extent, e.g. cluster.lower = std::min(cluster.lower, min(gamma_energies)) and cluster.upper = std::max(cluster.upper, max(gamma_energies)) (subject to lowest/highest_energy and a small margin), so no gamma counted in the significance test can fall outside the emitted ROI. Optionally, re-evaluate counts_in_region/significance over the final [cluster.lower, cluster.upper] before keeping the ROI so the keep decision and the emitted window agree. Treat as low priority/cleanup since shipped defaults do not trigger it.

### [LOW] All rel-eff curves' gammas are pooled before clustering, so source and NORM lines can be merged into one ROI and the per-curve attribution / amplitudes are mixed
- **Where:** `src/FitPeaksForNuclides.cpp:3633-3680, 3771-3804, 4122-4131`  _(R5:clustering)_
- **Issue:** gammas_by_counts is a single flat list across all rel_eff_index (3633-3680). ClusteredGammaInfo only stores gamma_energies and gamma_amplitudes (header struct, FitPeaksForNuclides.cpp:199-205) with no record of which rel-eff curve / which source produced each line. In the multi-curve auto path (7012, source curve + NORM curve), a source gamma and a NORM gamma within cluster_num_sigma*sigma get clustered and merged together (4122-4131), and the dominant-gamma logic (3945-3953, max_element over mixed amplitudes) can pick a NORM line as the cluster's dominant peak, driving weighted_mean, effective_fwhm, step-continuum reference (should_use_step_continuum uses cluster.gamma_amplitudes max, 3476-3478) and merge decisions. There is no way downstream to know a ROI mixes curves.
- **Action:** Tag each pooled gamma with its rel_eff_index/source in gammas_by_counts and ClusteredGammaInfo (e.g. carry a parallel curve-index vector), then (a) record per-curve composition on each ClusteredGammaInfo so downstream code and debugging can see when an ROI mixes source+NORM, and (b) make the ROI-geometry statistics robust to the asymmetry: for weighted_mean, effective_fwhm, the step-continuum reference gamma, and the merge dominant-gamma, consider restricting to source-curve gammas (or computing source-only and combined variants) so that turning NORM fitting on/off does not change source ROI placement/width or the step-vs-linear continuum choice. At minimum, verify and add a regression test that ROI bounds and continuum type for a fixed source spectrum are stable with NORM fitting toggled on vs off. This is low-cost and directly supports the goal of fitting reliably whether or not a background is present.

### [LOW] Greedy ordering by expected counts makes clustering, merge and break results depend on the rel-eff estimate quality and the source-fit order
- **Where:** `src/FitPeaksForNuclides.cpp:3692-3699, 3717-3733`  _(R5:clustering)_
- **Issue:** Clusters are seeded in descending expected-counts order; whichever gamma is processed first claims its +/-num_sigma window and erases neighbors from gammas_by_energy (3781), so a weaker neighbor can never seed its own (possibly better-centered) cluster. Expected counts = numPerSecond*rel_eff (3677). When the rel-eff estimate is poor (early manual stage, or extrapolated at energy extremes where manual_rel_eff is clamped, 5627-5635), the ranking — and thus which line becomes a cluster center vs an absorbed satellite — can change, shifting ROI centers. Tie-break is energy (3698) so it is deterministic, but it is sensitive to small rel_eff perturbations.
- **Action:** Treat as a low-priority robustness improvement, not a bug fix. Options, in rough order of value: (1) Decouple cluster window centering from the seed: after greedily grouping, re-derive each cluster's window from its amplitude-weighted-mean and effective FWHM (already computed at 3862/3899) in a fixed-point / second pass, so the final bounds do not depend on which member seeded the group. (2) Make seeding insensitive to rel_eff at the energy extremes by either not clamping or by tie/near-tie-breaking on energy when expected counts are within a small relative tolerance (the count comparison at 3696 currently uses exact inequality). (3) For diagnostics, log when a re-clustering pass changes a line's center/satellite role between the manual and auto stages so the GA fitness can penalize ROI instability. No change is required for correctness.

### [LOW] Initial significance numerator (amplitude*0.7607) vs denominator (full-FWHM data) are inconsistent windows
- **Where:** `src/FitPeaksForNuclides.cpp:2803-2813`  _(R6:observable-peaks)_
- **Issue:** peak_contrib = amplitude*0.7607 approximates Gaussian area within +/-1 FWHM, but data_area = foreground->gamma_integral(mean-fwhm, mean+fwhm) is the TOTAL counts (signal+continuum+background) over +/-1 FWHM. So significance = signal_in_window / sqrt(total_in_window). For a peak on a large continuum this is fine as S/sqrt(B), but `amplitude` here is the peak's full area (the Gaussian normalization), and 0.7607 is applied to it; if `amplitude()` already returns peakArea (full area) this is the +/-1FWHM fraction, consistent. The fragility is that data_area is clamped to max(...,1.0) (2813): for a near-empty high-energy ROI (data_area<1), sqrt floors at 1.0, so significance == peak_contrib, which can pass 2.25 on a few-count peak, admitting spurious high-energy peaks where Poisson noise is tiny but counts are negligible.
- **Action:** Treat as a minor robustness nit, not a correctness bug; do not retitle it as a window inconsistency. If touched, make the initial pre-filter floor smoother to avoid a step in GA's response surface near data_area~1, e.g. add a small expected-background/continuum term in the denominator (sqrt(data_area + eps) or sqrt(max(data_area, expected_continuum_in_window))) instead of a hard clamp to 1.0, so a sub-1-count window yields a near-zero (not floored) significance. Lower priority than the user's stated themes (knee, physical-earlier, init-rel-eff, fwhm-constraint). Optionally add a developer-check assert/log when a peak passes the initial filter only because of the 1.0 floor (data_area < 1.0) and later fails the final filter, to quantify how often this leaks via refit-failure (lines 2973-2974).

### [LOW] compute_significance_in_region: synthetic signal compared to sqrt(data) that includes the same signal
- **Where:** `src/FitPeaksForNuclides.cpp:3320-3335`  _(R6:observable-peaks)_
- **Issue:** significance = sum_synthetic / sqrt(sum_data). sum_data is raw data over the check window, which contains the real peak counts too, not just background/continuum. The proper S/sqrt(B) would use background-only in the denominator. Using total data inflates the denominator for strong peaks, so the metric saturates: a very strong peak and a moderate peak can yield similar 'significance' because both numerator and (sqrt of) denominator grow with amplitude. This is used both for breakpoint selection (find_synthetic_minima, has_significant_peak_between) and indirectly to decide whether to break a wide cluster.
- **Action:** Treat as a heuristic-robustness improvement, not a hard bug. (1) Make the noise estimate background/continuum-aware instead of raw total counts: estimate B in the check window from a local linear continuum (or subtract a live-time-scaled background estimate) and use S/sqrt(max(B, small_floor)); this removes both the strong-peak saturation and the with/without-background asymmetry. (2) For find_synthetic_minima, consider that the ranking key is computed at valleys, so the higher-value fix there is ensuring the valley significance reflects the real dip depth in data, not just total counts. (3) If keeping the current form, document that this stage is intentionally foreground-only and that break_peak_significance_threshold (2.0) is calibrated against S/sqrt(S+B), so the GA tuning per CoarseResolutionType already partly compensates. Low priority relative to the four named themes (knee, physical-earlier, init-rel-eff, fwhm-constraint).

### [LOW] NORM fit forces background subtraction off, so NORM-curve and subtraction are mutually exclusive
- **Where:** `src/FitPeaksForNuclides.cpp:8494-8519`  _(R7:roi-norm-escape)_
- **Issue:** When FitNormBkgrndPeaks/DontUse is set, long_background=nullptr (8497) so no background subtraction occurs, and a second FramPhysicalModel NORM rel-eff curve is added (6624-6769). The GA selects this purely via BackgroundMode (PeakFitImprove.cpp:1748; NuclideConfig_GA.cpp:301), independent of whether src_info.long_background exists. Thus there is no path that BOTH subtracts a measured background AND models residual NORM lines; and the NORM source list is the fixed 5 nuclides regardless of what the actual background spectrum contains (e.g. Cs137/Co60 ambient lines, or a NORM line absent from this site).
- **Action:** Treat as a design-limitation / improvement opportunity, not a bug. Two concrete options, both relevant to the stated themes: (1) Allow a combined path that BOTH subtracts a measured long_background AND keeps the NORM rel-eff curve for residual NORM lines (remove the unconditional long_background=nullptr at 8497 when a real background is present, and let the second FramPhysicalModel curve model whatever subtraction leaves behind). (2) Make the NORM source set adaptive: rather than the fixed five chains in get_norm_sources (761-787), seed candidate NORM/ambient nuclides from peaks actually present in the supplied background spectrum (e.g. Cs137 661, Co60 1173/1332) so background-present spectra do not silently omit real ambient lines. If left as-is, document at 8497 that subtraction and NORM-modeling are intentionally exclusive and that ambient non-NORM-chain lines are not modeled in the NORM-fit path.

### [LOW] is_near_strong_norm_gamma table documented as sorted but is not; energies also drift from current data
- **Where:** `src/FitPeaksForNuclides.cpp:863-909`  _(R7:roi-norm-escape)_
- **Issue:** The comment says 'Keep sorted ... helps future binary-search optimizations' but the array is not sorted (e.g. 84.94 then 63.29; 277.36 before 241.99). A future binary_search refactor would silently miss lines. Linear scan is correct today, so impact is latent. Separately, a few entries are hand-curated approximations (e.g. 129.07 'Ac-228/Pb-214 mix', 270.24) that may not match the DB energies used elsewhere.
- **Action:** Comment-only fix (no behavior change): remove or correct the "Keep sorted ... future binary-search optimizations" comment at src/FitPeaksForNuclides.cpp:861-862 to accurately describe the array as grouped-by-nuclide and consumed by a linear scan, OR genuinely sort the array ascending by energy so the comment becomes true and a binary-search refactor would be safe. Do not treat the hand-curated/DB-drift concern as a defect. Given the linear scan is correct, this is low priority.

### [LOW] assign_escape_peak_relationships re-runs on three peak vectors that may share PeakDef state; significance uses peakAreaUncert that can be unset
- **Where:** `src/FitPeaksForNuclides.cpp:1607-1758, 8161-8163`  _(R7:roi-norm-escape)_
- **Issue:** Called separately on result.fit_peaks, result.observable_peaks, and result.uncombined_fit_peaks. Significance test uses peak_area/peakAreaUncert with the guard (uncert>0?...:0) (1713,1740); if uncert is 0/unset (can happen for combined peaks or peaks where uncertainty wasn't populated), significance is 0 and the escape relationship is silently not assigned, so fit_peaks vs uncombined_fit_peaks can get different escape attributions. Parent match tolerance is a flat 2.0 keV (1648,1708,1735) regardless of detector resolution.
- **Action:** Make the significance test robust to an unpopulated area uncertainty instead of treating it as significance 0. Options: (a) fall back to a Poisson estimate of the uncertainty (sqrt(area), optionally plus continuum) when peakAreaUncert() <= 0, so escape labeling still works when covariance failed; or (b) gate the whole significance test on a present uncertainty and otherwise use an absolute area threshold. Also fix combine_peaks (FitPeaksForNuclides.cpp:2457-2476) to skip or zero-out negative (sentinel) input uncertainties so it never manufactures a positive uncert from -1 values, which removes the fit_peaks-vs-uncombined_fit_peaks divergence. Consider scaling the 2.0 keV parent/escape match tolerance by the detector FWHM at the relevant energy. Severity is low: the path only triggers on HPGe with fit_norm_peaks, for the single Th232/Tl208 2614 keV candidate, and only when covariance computation fails or escape peaks get combined.

### [LOW] combine_peaks averages sigma/mean by area but adds continuum-correlated uncertainties in quadrature, and assumes shared continuum
- **Where:** `src/FitPeaksForNuclides.cpp:2434-2490`  _(R7:roi-norm-escape)_
- **Issue:** setPeakAreaUncert(sqrt(sum_uncert^2)) treats the merged sub-peaks' area uncertainties as independent, but peaks sharing a ROI/continuum are strongly correlated, so the combined uncertainty is underestimated. combined sigma is area-weighted average of sub-peak sigmas, which understates the width of two separated peaks merged into one. The assert(peak->continuum()==cont) (2452) will fire (or throw in release via the following if) if group_peaks_by_roi ever passes peaks with differing continuum pointers; combine_overlapping_peaks_in_rois groups by continuum so this holds today, but is fragile.
- **Action:** If desired (low priority), fix combine_peaks (src/FitPeaksForNuclides.cpp:2474-2486) to (a) compute the merged sigma including the between-means spread, e.g. sigma_combined = sqrt( sum(area_i*(sigma_i^2 + (mean_i - mean_combined)^2)) / total_area ), and (b) for the merged area uncertainty, account for positive correlation among co-continuum peaks rather than pure quadrature (a conservative bound is the linear sum sum|uncert_i|, or fold in the fit covariance if available). The assert/throw pair (2452-2454) can stay as-is; it is a correct invariant guarantee, not a bug. No change is needed for GA correctness, because the GA scores result.observable_peaks (NuclideConfig_GA.cpp:1131) from compute_observable_peaks, which never passes through combine_peaks; combine_peaks only affects result.fit_peaks (the display/reporting set, single call site at FitPeaksForNuclides.cpp:7546).

### [LOW] shrink_rois_for_interfering_peaks gathers all gamma energies globally, so a source gamma in one ROI suppresses shrinking in another
- **Where:** `src/FitPeaksForNuclides.cpp:433-488`  _(R7:roi-norm-escape)_
- **Issue:** all_gamma_energies is collected across ALL rois_and_gammas (433-439). An interfering auto-search peak is exempted from shrinking if it is within cluster_num_sigma of ANY source/NORM gamma anywhere in the spectrum, not just gammas in the ROI it overlaps. With the NORM second curve active, the NORM gamma list is large (many chain lines), so a genuine interfering peak that happens to sit near any NORM line is treated as 'ours' and never trimmed.
- **Action:** Optional low-priority hardening: in shrink_rois_for_interfering_peaks, scope the exemption test to the gammas of the ROI(s) the peak overlaps (or to gammas whose own surviving ROI overlaps the peak) instead of the global all_gamma_energies, so a NORM/source line in an unrelated region cannot suppress trimming of a different source ROI. If this global match is kept intentionally, add a comment noting it is deliberate and that the tolerance equals the clustering half-width so the cross-ROI leakage is bounded. Do NOT tie any of this to background presence; the relevant switch is the FitNormBkgrndPeaks option, and the same scoping consideration applies to the merge-prevention code at lines 3920-3985 for consistency.

### [LOW] GAOptimized chi2-cap gene can produce FramPhysicalModel with order>0, tripping a debug assert; release self-corrects silently
- **Where:** `target/peak_fit_improve/NuclideConfig_GA.cpp; src/FitPeaksForNuclides.cpp:NuclideConfig_GA.cpp:540-542,645-648,814-815; FitPeaksForNuclides.cpp:5773-5777`  _(R8a:GA-wrapper)_
- **Issue:** rel_eff_eqn_type gene ranges 0..4 (FramPhysicalModel=4) and rel_eff_eqn_order gene ranges 0..6, independently. genes_to_settings clamps each but does not force order=0 when type==FramPhysicalModel. The engine asserts (FitPeaksForNuclides.cpp:5774) that order==0 for FramPhysicalModel before self-correcting to 0 at line 5777. In a debug/asserts build the GA will abort the first time it samples that combination; in release it silently forces order=0.
- **Action:** Couple the two genes so FramPhysicalModel always maps to order 0 (and ideally is excluded from order mutation/crossover). Cleanest fix is in NuclideConfig_GA.cpp genes_to_settings (around 540-542): after clamping, if config.rel_eff_eqn_type==FramPhysicalModel set config.rel_eff_eqn_order=0. This both prevents the debug abort and collapses the redundant genome region at the source so the GA does not see spurious diversity. Optionally also normalize p.rel_eff_eqn_order=0 in init_genes/mutate/crossover when type==Fram to keep stored genomes canonical. Leave the engine assert/self-correct at FitPeaksForNuclides.cpp:5774-5777 as a defensive guard for non-GA callers.

### [LOW] Best-individual re-evaluation for fg/bg breakdown is non-deterministic and can disagree with stored best cost
- **Where:** `target/peak_fit_improve/NuclideConfig_GA.cpp:994-1010`  _(R8a:GA-wrapper)_
- **Issue:** When sm_do_background_fit_trial is on, SO_report_generation re-runs ns_ga_eval_fcn(best_config) to recover fg/bg components (line 1000). The fit itself is deterministic for a fixed config, but the parallel summation path (PeakFitImprove.cpp:1809-1825) sums in spectrum-index order so that is stable; however the reported best_fg/best_bg come from this RE-EVALUATION, not from the cached middle_costs of the actual best individual. If genes_to_settings(best_genes) is byte-identical this matches, but any nondeterminism inside fit_peaks_for_nuclides (e.g. Ceres threading) would make the reported breakdown inconsistent with last_generation.best_total_cost. Also costs one extra full-dataset eval per new-best generation.
- **Action:** Replace the re-evaluation block (NuclideConfig_GA.cpp:994-1008) by reading the cached breakdown directly from the best chromosome: when best_yet, set sm_best_fg_cost/sm_best_bg_cost from last_generation.chromosomes[last_generation.best_chromosome_index].middle_costs.objective_fg/objective_bg, and fix/remove the misleading comment at line 988. Guard with a check that best_chromosome_index >= 0 and < chromosomes.size(). This makes the reported breakdown reconcile exactly with best_total_cost and removes one redundant full-dataset eval per improvement. Low priority since it only affects diagnostic logging in the optimization harness, not peak-fit correctness or GA selection.


---

# 2026-07-19 — Phase 3 · R1 step 2 (cont.): global SNIP continuum wired into the remaining gates

Working tree, NOT committed.  Base commit `49624ea2` (clean tracked tree; concurrent session's
determinism fix already landed).  Consumer #1 (keep/drop gate) and the provider
(`detail::GlobalContinuumEstimate` / `make_global_continuum`, built once from `orig_foreground` at
`fit_peaks_for_nuclide_relactauto` ~:7372) were already landed.  This change finishes the per-consumer
threading the tingly-abelson plan §"Step 3 — R1 step 2" specifies, EXCEPT consumer #5 (which stays local
by design).  Each consumer landed separately with the full `test_fitPeaksForSources` suite green between.

Pattern (all three): prefer `global_continuum->integral(...)` when `global_continuum && ->valid()`, else
fall back to the EXACT prior local→gross estimate.  Behaviour is therefore byte-identical wherever the
provider fails to build (`valid()==false`), so this cannot regress the invalid-provider path.

- **Consumer #3 — edge-ROI restore fallback** (`src/FitPeaksForNuclides.cpp` ~:8726, the initial-edge-ROI
  restore inside the refinement re-cluster).  The weakest local estimator (no signal subtraction).  The
  `b_est` for the restore detection-z now comes from `global_cont` when valid, else the former
  `estimate_local_continuum`→`gamma_integral` chain.  Only affects the branch where the auto-search peak
  lacks an area uncertainty.  `global_cont` in scope directly (same fn).
- **Consumer #2 — step-gate Gate-1 dominance** (`cluster_gammas_to_rois` ~:5498).  ONLY the Gate-1
  dominance `b_est = step_cont.integral(win_lo,win_hi)` now prefers `settings.global_continuum`.  The
  `if(step_cont.valid)` guard and Gate-2 `step_cont.sideband_asymmetry_z()` stay LOCAL (they need raw
  sideband counts SNIP does not provide); the "no gross fallback for a step" invariant is preserved
  (fallback is the local `step_cont.integral`, still gated on `step_cont.valid`).
- **Consumer #4 — `find_clean_gap_between` continuum anchor** (`detail::find_clean_gap_between`, decl
  `InterSpec/FitPeaksForNuclides.h` ~:173 + fwd-decl of `GlobalContinuumEstimate`; def `src/...` ~:855).
  New defaulted trailing `const GlobalContinuumEstimate *global_continuum = nullptr`; the window
  continuum `c_est` prefers it when valid (predicted-signal tail term unchanged).  Wired at BOTH callers:
  the re-cluster path (`cluster_gammas_to_rois` ~:4978 → `settings.global_continuum`) and, via a
  defaulted param added to file-local `merge_rois`, the in-solve initial merge (`merge_rois` call at
  ~:8109 → `global_cont`).  The two pre-solve `merge_rois` sites (`estimate_initial_rois_without_peaks`
  ~:6085 and `fit_peaks_for_nuclides` ~:10396) keep the default nullptr — `global_cont` is not built
  until inside `fit_peaks_for_nuclide_relactauto`, so no provider exists there; using `orig_foreground`
  as the canonical B(E) basis is consistent with the keep-gate wiring and the plan's "provider stays on
  original cal" decision.

**Verification:** baseline (`49624ea2`) full suite GREEN (all 7 suites); then #3, #2, #4 each built and
full `test_fitPeaksForSources` GREEN in turn (SmokeTests / IdempotencyTests / OptionBehaviorTests /
TrinititeSequential / AdditionalSpectra / BystanderDegradation / StatisticalDetailHelpers — "*** No
errors detected").  Trinitite (NORM-heavy: exercises keep/step/merge on K40/Eu152 ROIs) is the sensitive
suite and stayed green across all three.  Adversarial re-review + optional eval/gallery A/B still to run.

**GA note — re-tune RECOMMENDED (not gene-value change):** R1-step2 changes the *source* of B(E) for the
keep-gate, step-gate Gate-1, edge-ROI-restore and clean-gap-merge decisions (global SNIP vs local
two-sideband).  No gene semantics/values were edited, but the keep/step/merge thresholds now see a
different B, so a per-det-type GA RE-EVALUATION is warranted (deferred, big machine); the manual eval
loop is the interim calibration.  Eval determinism is fixed (`49624ea2`), so a same-batch A/B is now
trustworthy when run.

## R1-step2 (#2/#3/#4) adversarial review — 2026-07-19: NO CONFIRMED BUGS
Independent subagent trace of the +43/-13 diff. Verdict: no bugs introduced; safe to compile+run. Verified:
merge_rois threading correct (`detail::` qualification; all 3 call sites compile — two default to nullptr,
:8109 passes `global_cont.valid()?&global_cont:nullptr` on the SAME orig_foreground basis, consistent);
fallback byte-equivalent for all three consumers on the invalid-provider path (Gate-2 asymmetry + step_cont.valid
guard intact for #2; only the peakAreaUncert<=0 branch touched for #3; c_est chain prepend-only for #4);
every deref guarded by `&& ->valid()` / `.valid()`; `global_cont` lifetime sound; header forward-decl clean
(incomplete type OK for the pointer param, default arg declared once); integral() contract matches the prior
LocalContinuumEstimate::integral.
- Carried note (PLAUSIBLE, low, PRE-EXISTING, by design — NOT introduced here): `global_cont` is built from
  `orig_foreground` (original cal) but the working `foreground` is cal-advanced inside the refinement loop
  (~:8562 when fit_energy_cal, default true), so the SNIP integral is queried on original-cal binning while
  the local fallbacks + data use the fitted cal. This ALREADY affects the landed keep-gate (consumer #1) and
  is the plan's explicit "provider stays on original cal; sub-keV gating error; do not recompute" decision;
  R1-step2 merely extends the same approximation to the step/edge/merge gates. SNIP is smooth => first-order
  cal-insensitive; low impact. Optional future refinement (NOT an action item): rebuild global_cont on cal advance.
- Trivial (not a bug): find_clean_gap_between still computes local `cont` even when the global is used (wasted work only).

## R1-step2 cal-frame fix + accuracy A/B — 2026-07-19
- **Cal-frame fix (LANDED, working tree):** `global_cont` is now re-stamped to the fitted energy cal whenever
  the refinement loop advances the cal (`fit_peaks_for_nuclide_relactauto`, at the `foreground = new_foreground`
  site). The fit only RE-LABELS the cal (copies channel counts, `set_energy_calibration`), so the SNIP's
  per-channel counts are invariant — we copy `global_cont.snip` + `set_energy_calibration(fitted_cal)` and point
  `.foreground` at the advanced working foreground (cheap; no SNIP recompute). Closes the review's F1 cal-basis
  note for ALL consumers (incl. the already-landed keep-gate #1). `global_cont` made non-const for this. Full
  `test_fitPeaksForSources` GREEN.
- **Adversarial review of the R1-step2 diff: NO CONFIRMED BUGS** (recorded above).
- **Accuracy A/B (R1-step2 ON vs OFF, `--config-genes default`, 300s, holdout-frac 0.03, seed 20260703, 8 threads;
  determinism fixed at 49624ea2 so deltas are real):**
  | det | OFF avg cost | ON avg cost | delta | n | mech-fail | empties |
  |-----|-------------:|------------:|------:|--:|----------:|--------:|
  | HPGe | -8.50692 | -8.50704 | -0.0001 (neutral) | 127 | 0 | 1/1 |
  | NaI  |  3.91129 |  3.93113 | +0.020 (+0.5% worse) | 292 | 0 | 0/0 |
  Interpretation: R1-step2 is accuracy-NEUTRAL (HPGe) to marginally-negative (NaI +0.5%) on the current
  UNVALIDATED default config, with ZERO reliability impact. The gates' thresholds were tuned for the LOCAL B(E)
  and now see the GLOBAL B(E) → the accuracy verdict is gated on the DEFERRED per-det-type GA re-tune. Value is
  architectural (single shared B(E); prerequisite for R2/R4 and a clean re-tune). Recommendation: KEEP.
  (Logs: run_r1_eval/{HPGe,NaI}_r1{off,on}.log. NaI r1-on was re-run standalone after the first attempt was
  killed mid-run when its launching subagent was reaped — lesson: run long evals as a harness-owned job.)


---

# 2026-07-20 — R6 debt closure, R2 bounded rescue, and R4 shadow evidence

Working tree, **NOT committed**, based on `47ee4a5c` (R1 step 2 is already committed at that HEAD).
No GA tuning or production ambient/background/floating-interferer search was enabled.

## R6 debt closure — LANDED

- Production discovery remains limited to foreground-confirmed strong NORM lines.  Cs137/Co60 ambient
  sweeping, dedicated background discovery, and unattributed floating interferers remain disabled.
- Candidate selection is deterministic and bounded to at most two **usable** nuisance nuclides per call,
  ranked by maximum confirming significance and stable nuclide identity.  Each selected parent must
  contribute a fitted confirming line at `z >= 5`; peak-object presence alone is insufficient.
- A successful source-only incumbent is preserved before nuisance augmentation.  The challenger is accepted
  only after a common-domain solve, requested-anchor preservation outside affected neighborhoods, and a
  significant affected-ROI Poisson-deviance likelihood improvement at `config.roi_significance_z`; otherwise
  the incumbent is retained with a warning.
- Automatic nuisance siblings are excluded from ordinary source re-clustering.  Accepted nuisance lines are
  kept in the private model through ROI combining and observable LM refitting, cannot combine/cull across
  source-vs-nuisance provenance, and are removed only from the three final public peak vectors.
- `PeakFitLM::fit_peaks_in_spectrum_LM` gained a defaulted close-pair-removal predicate; default behavior is
  unchanged for every other caller.  `compute_observable_peaks` now distinguishes peaks that **must reach**
  the honest LM refit from nuisances that **must remain** through final private-model filtering.
- Observable re-absorption question: code trace confirmed that a source-blind significance removal could drop
  a nuisance before the next observable iteration, permitting a source-only refit to absorb its counts.  The
  filter-late path above closes that mechanism.  The strong K40/weak Eu152-1457 explicit-joint regression now
  stays within the larger of 3 combined sigma or 20%, while K40 is absent from all public vectors.
- Tests include the raw two-parent R6 transaction, bystander suppression, explicit joint comparison, and a
  six-fit `{Eu152,Cs137}` / reversed-order / `{K40,Eu152,Cs137}` stress matrix in raw and supplied-background
  modes.  No zero-peak collapse and requested anchor areas remain within 3 sigma/20%.

**Final independent R6 adversarial review:** no confirmed functional bug.  The reviewer traced discovery,
ranking/cap, curve/ROI construction, source-only incumbent, augmented admission, refinement, observable
provenance, and final filtering after all fixes.  Deferred recommendation: keep ambient/background/floating
discovery parked until a separately conditioned representation is designed and stress-tested.

## R2 bounded fit-then-prune rescue — LANDED

- Clustering now uses provenance-preserving predicted-gamma records (energy, expected counts, source, curve
  index).  Normal accepted ROI semantics remain unchanged.  Rejected clusters are recorded as marginal only
  when they pass the existing count floor, have significance in `[0.7*keep_z, keep_z)`, and fail the normal
  keep threshold.
- One pass runs after the successful incumbent and guarded R6 transaction.  Candidates are re-predicted from
  fitted activities/rel-eff curves, rebuilt by the same extent/continuum-order clustering machinery, ranked
  deterministically, and limited to four actual provisional ROI Wilks trials.
- Guards exclude accepted/rejected/previously-filtered overlap, rel-eff extrapolation (including fewer than two
  distinct fitted source energies), unfit auto peaks, every detected/warned R6 line, and strong unmodeled NORM
  lines within one FWHM.  Rejected ranges are never reconsidered.
- Passing candidates enter one refinement solve.  The entire challenger—option preparation, solve, anchors,
  and post-solve evaluation—is locally transactional; any returned failure or exception retains the incumbent.
  A forced post-solve evaluation exception proves this rollback path.
- Rescued source peaks bypass only the coarse observable prefilter/internal LM significance cull so the earlier
  ROI Wilks decision reaches an honest refit.  They still must pass the ordinary final post-refit source
  significance threshold.  The final insignificant-ROI filter remains authoritative.
- Fixed rails remain file-local constants: rescue fraction `0.7`, guard `1.0 FWHM`, maximum four ROIs.
- Tests cover unchanged marginal classification, causal enabled/disabled recovery, post-solve exception
  rollback, strong-interferer non-overfit, absent Eu154/no resurrection, source-order/raw-vs-background fitted
  area invariance, and per-mode observable order invariance.

**Same-batch R2 ON/OFF (`default`, 300 s, 3% eval holdout, seed 20260703, 8 threads):**

| Detector | OFF avg cost | ON avg cost | Missed OFF/ON | Spurious OFF/ON | Mechanical failures | Legitimate empties | Stable-ID rows changed |
|---|---:|---:|---:|---:|---:|---:|---:|
| HPGe (127) | -8.62109 | -8.62109 | 37 / 37 | 11 / 11 | 0 / 0 | 1 / 1 | 0 / 127 |
| NaI (292) | 3.83685 | 3.83685 | 77 / 77 | 118 / 118 | 0 / 0 | 0 / 0 | 0 / 292 |

Stable-ID comparison found zero differences in cost, status, miss fraction, observable count, fit count, or
error.  Interpretation: rescue is intentionally dormant under current default thresholds on this holdout;
the causal synthetic test proves it activates when a real line is deliberately moved into the marginal band.
No reliability regression.  **Per-detector GA re-tuning/re-evaluation is warranted after gate semantics freeze,
but remains deferred.**

**Final independent R2 adversarial review:** no confirmed bug after re-reviewing normal-clustering semantics,
support spans, rejected ranges, four-trial cap, whole-challenger exception safety, causal controls, interferer
guards, and the fitted-vs-observable significance contract.

## R4 SNIP-guided joint boundaries — SHADOW PROTOTYPE ONLY; DO NOT ENABLE

- Shadow mode is gated by `PEAKFIT_ROI_SHADOW_TSV` / `INTERSPEC_ROI_BOUNDARY_SHADOW` and never replaces
  production ROIs.  Initial manual/fallback calls build a separate shadow-only SNIP estimate; invalid SNIP or
  infeasible partitions fall back diagnostically while production remains legacy.
- Candidate intervals use atomic accepted gamma groups, source cores, legacy bounds, gap midpoints,
  SNIP valleys/curvature, and unfit-peak exclusion edges.  Linear, quadratic, flat-step, and linear-step
  continua are fit to SNIP with foreground Poisson variance; the existing maximum FWHM width is only a
  catastrophic backstop.  Diagnostics include stable file-based spectrum IDs, old/proposed geometry, channels,
  FWHM width, continuum family/mismatch, reasons/conflicts, and foreground/SNIP/fitted-continuum profiles.
- Evaluated separately on HPGe (127 spectra), NaI (292), LaBr (105), and CZT (139), all with zero mechanical
  failures.  Reproducible TSVs and class/combined HTML review packages are in
  `target/peak_fit_improve/r4_shadow_evidence/`.

| Class | Diagnostics | Valid | Fallback | Changed bounds | Actual splits | Actual merges | Multi-gamma-group intervals | Max proposed width |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HPGe | 312 | 28 | 284 | 23 | 0 | 0 | 20 | 14.98 FWHM |
| NaI | 665 | 234 | 431 | 218 | 1 | 0 | 232 | 19.94 FWHM |
| LaBr | 247 | 80 | 167 | 70 | 0 | 0 | 78 | 19.82 FWHM |
| CZT | 312 | 81 | 231 | 68 | 0 | 0 | 77 | 19.71 FWHM |

**Verdict:** useful negative evidence, not a production candidate.  The additive normalized interval score
charges average mismatch once per ROI and therefore has a structural bias against splits; only one actual
split and no actual partition merges were proposed across all detector classes.  The diagnostic reason named
`merge` currently means that one interval covers multiple atomic modeled-gamma groups, which usually already
shared a legacy ROI; the multi-group counts above must not be read as legacy ROI merges.  Most fallbacks were
`no feasible non-overlapping SNIP partition`.  The next revision
must score explicit ROI-vs-clean-gap decisions on a common comparison domain with a calibrated ROI-complexity
penalty, fix multi-legacy-ROI visualization/lineage, and add the missing raw-data valley/normalized slope
candidates.  Do not interpret current score deltas as likelihood/AIC improvements and do not promote any R4
threshold to a gene.  Production enablement remains a separate user-review decision.

**Independent R4 adversarial review:** production isolation and the major initial-visibility, atomic-lineage,
unmodeled-conflict, stable-ID, and visual-profile issues were confirmed fixed.  The final evidence audit
verified every diagnostic/valid/fallback/changed/max-width count, caught and corrected the multi-group-vs-
actual-merge label above, and supported the conservative `DO NOT ENABLE` verdict.

## Final verification

- `cmake --build target/peak_fit_improve/build_vscode`: GREEN after final source changes.
- `cmake --build target/testing/build_test --target test_fitPeaksForSources`: GREEN.
- Final full `test_fitPeaksForSources`: **33 test cases, GREEN (`*** No errors detected`)**.
- Focused StatisticalDetailHelpers (8), SmokeTests (7), TrinititeSequential (7), BystanderDegradation (2),
  raw R6 transaction, R2 causal/rollback/invariance/interferer tests: GREEN during development.
- `git diff --check`: clean.  Nothing staged or committed.

---

# 2026-07-20 — caller-controlled R6 isolation and manual-tuning strategy

- Added `FitSrcPeaksOptions::DisableAutoInterfererFit` as a backward-compatible opt-out. The default
  still permits guarded R6 nuisance fitting; the flag skips the complete automatic detection,
  warning, nuisance-curve, and nuisance-ROI path. Explicit requested sources and explicit NORM-fit
  options are unchanged. This is intentionally independent of whether a background pointer exists.
- Added `PeakFitImprove --disable-auto-interferer-fit`. Foreground scoring, optional background-
  false-positive trials, checkpoint compatibility metadata, and generated fit HTML all use and label
  the same setting, so a source-only tuning run cannot silently render R6-enabled plots.
- Extended the raw R6 transaction regression: default options retain the strong K40 nuisance privately,
  while the opt-out fit contains only the two requested parents and still returns their anchors.
- Planning decision: no immediate GA tune. First manually tune low-, medium-, and high-resolution
  classes separately with R6 disabled, using scalar components plus per-spectrum visual review of ROI
  geometry, continuum family/plausibility, step use, peak grouping, areas/means, misses/spurious peaks,
  and stability. Generate an all-spectrum HTML gallery for each stable candidate, obtain user review,
  revise the implementation/metric/parameters, and iterate. Re-enable and assess R6 as a distinct
  later layer. A GA is reserved for distant optional fine-tuning after the pipeline and judgment
  metric are demonstrably reliable.
- Verification: engine and test targets rebuilt cleanly; the focused raw R6 enabled/disabled
  transaction passed; `PeakFitImprove --help` exposes the new switch; and the authoritative full
  `test_fitPeaksForSources` run completed all 33 cases with `*** No errors detected`.

---

# 2026-07-20 — manual-review HTML reporter audit and rebuild

## Audit verdict

The former `NuclideConfigEval` HTML path was not trustworthy enough for manual tuning.  Command-line
scoring and HTML generation independently refit every spectrum, so the displayed fit was not guaranteed
to be the result that produced the printed objective and report generation doubled the expensive fit work.
Non-success results were assigned obsolete report-local accounting, then omitted from the gallery; empty
and failed cases could therefore disappear.  The page plotted public observable peaks without identifying
that representation, did not expose the combined/uncombined/private representations, had weak supplied-
background labels, and did not show final ROI geometry, continuum choice, exact objective matching, or
component costs.  Filesystem enumeration and completion-order floating-point reduction also made selection
or aggregate output sensitive to ordering/thread count.

The data path is now explicit and single-pass:

1. `NuclideConfigEval` deterministically sorts selected `PrecomputedNuclideData` records.
2. Each worker writes one index-addressed `SpectrumEvaluation`, including the complete `PeakFitResult`,
   status, warnings/errors, exact shared foreground score breakdown, and optional background trial.
3. Component totals are reduced in stable spectrum order after workers join.
4. The HTML writer consumes that exact cache and refuses missing, duplicate, or mismatched records.

The scalar objective and report use `PeakFitResult::observable_peaks`.  Candidate matching is shared with
the objective: a found mean inside the truth ROI, or its one-sigma interval overlapping that ROI; escape-peak
correction and the first unexpected 508--514 keV exemption are preserved.  FinalFit's separate 1.5-sigma
nearest-truth area comparison remains unchanged.  The gallery's “largest area pull” is consequently labeled
report-only rather than being presented as an objective term.

The production model overlay is constructed from public combined `fit_peaks`; its continua and bounds are
the final shared `PeakContinuum` objects on those peaks.  Objective-observable and combined-public means are
marked separately.  Counts for `uncombined_fit_peaks`, `observable_peaks`, private
`solution.m_peaks_without_back_sub`, and final solver ROI ranges are disclosed.  Private solver means are not
plotted because they are expressed in the fitted-calibration frame and could be displaced over the original
foreground calibration.  R4 shadow output is never rendered as production.

Supplied background is displayed as the source foreground plus background scaled by
`foreground_live_time/background_live_time`, with an optional foreground-minus-scaled-background curve and
explicit live-time/scale labels.  The fitted production total remains the fit's own public model; the report
does not imply that a report-side net curve was fitted.  `DisableAutoInterfererFit` is applied once to the
cached scoring fit and is prominently recorded in the header, eliminating scoring/rendering disagreement.

## Reporter implementation

- Every selected spectrum receives a stable path-qualified FNV-1a ID and anchor, in deterministic canonical
  order.  Legal repeated metadata in different corpus paths is disambiguated; exact duplicate inputs are
  rejected as an invariant violation rather than silently overwriting an entry.
- Successes, legitimate successful empties, and every non-Success mechanical failure are counted separately.
  Failed partial results are visible but explicitly called diagnostic, with complete status/error/warnings.
- Report metadata records detector class/models, exact filters, genes, background/R6 mode, thread count,
  invocation, static/data directories, build stamp, timestamp, objective totals, major components, and
  missed/spurious counts.
- A sortable/filterable index links all spectra and ranks scalar cost, report-only area pull, misses,
  spurious peaks, ROI FWHM/channel width, peaks per ROI, step use, fitted-continuum/SNIP mismatch, and
  failure/empty status.  The page states that diagnostic rankings are not necessarily objective terms.
- Spectrum sections are collapsed and D3 charts are instantiated lazily.  Each includes metadata, complete
  errors, background scaling, measured/background/net curves, production fitted total/continuum, final ROI
  boundaries, truth/miss markers, combined and objective peak markers, ROI geometry/continuum diagnostics,
  fitted/public and exact-objective peak tables, truth matching, explicit spurious entries, background-trial
  details, and representation disclosure.  Hash navigation expands and scrolls to the requested section.
- User remarks can be attached to summary or detailed rows (including selected text), persist per report path
  in `localStorage`, and can be copied as an aggregated context-rich text block for the next Codex iteration.
- All untrusted text passes a shared HTML escaper; chart titles no longer interpolate source text into
  JavaScript.  The notes UI uses Safari-13-compatible DIV modals and avoids unsupported flex-gap/inset CSS.
- The old duplicate-refit reporter remains in the source only as disabled legacy code pending a later cleanup;
  the production path returns after the cached writer.

## Verification and artifacts

- Focused `PeakFitImproveReporterSelfTest`: GREEN.  It covers deterministic/path-disambiguated IDs, duplicate
  rejection, HTML escaping, ROI channel/FWHM measures, success/empty/failure coverage, and exact score-helper
  agreement fixtures.
- `cmake --build target/peak_fit_improve/build_vscode`: GREEN after the final source change (only existing
  OpenGA deprecation and SDK search-path warnings).
- `cmake --build target/testing/build_test --target test_fitPeaksForSources`: GREEN.
- Authoritative full `test_fitPeaksForSources`: **33 cases, GREEN (`*** No errors detected`)**.
- Final Detective-EX supplied-background/source-only smoke: 33 stable entries, 30 ordinary successes,
  2 legitimate empties, 1 authentic `FailedToSetupProblem` mechanical failure, 18 definite misses,
  1 spurious peak, multi-ROI/multi-peak ROIs, and stepped continua.  Automatic interferer fitting is visibly
  disabled.  The gallery contains 37 syntactically valid inline scripts and 33 unique spectrum anchors.
- Independent 8-thread and 1-thread runs produced identical ordered anchor hashes, summary-row hashes,
  objective `-315.786`, component values, and outcome counts.
- Headless Chrome inspection at 1440 px confirmed the summary's layout, labels, filtering controls, and
  failure/empty visibility.  Chrome's direct `file:` fragment screenshot remained blank even though the
  generated hash handler and anchors were validated; ordinary summary rendering and JavaScript parsing were
  successful.  This is retained as an inspection-tool limitation, not evidence that the link is absent.
- Review gallery: `target/peak_fit_improve/reporter_smoke_20260720/reporter_smoke.html`.
  Determinism comparators: `reporter_smoke_threads8.html` and `reporter_smoke_threads1.html` in the same
  directory.  Desktop screenshot: `reporter_smoke_summary.png`.

## Independent adversarial review and remaining limitations

The reviewer specifically audited plot representation, objective divergence, missing failures, ordering,
HTML/JavaScript injection, supplied-background scaling, and R4/R6 semantics.  Confirmed findings were fixed:
Safari `<dialog>`/CSS incompatibility, raw source text in a chart-script context, ambiguous public-vs-
observable peak overlays, production labeling of failed partial objects, missing background-trial details and
weight, completion-order reduction, non-disambiguated IDs, forced-log net data, absent fitted-mean markers,
opaque note context, and the fitted-calibration private-mean overlay.  No accidental R4 or R6 behavior change
was found.

Remaining constraints are explicit: the HTML embeds spectrum arrays up front even though charts instantiate
lazily, so very large cross-detector custom reports may still be memory-heavy; the nominal workflow should
continue generating one detector model per file.  Continuum Pearson RMS and fitted-continuum/SNIP mismatch
are report-only diagnostics, not likelihood claims.  Mechanical failures remain independent of scalar cost
rather than introducing a new tuning penalty.  Exact duplicate records pointing to the same corpus path are
treated as an input error.  No parameter tuning, GA run, R4 enablement, or interferer expansion occurred in
this increment.

## User-review follow-up: encoding and invisible production peaks

The first user inspection found two presentation defects that the initial headless summary screenshot did not
expose.  `D3SpectrumExport::write_html_page_header` did not declare UTF-8, so browsers could render em dashes
as `â€”`.  More importantly, the cached reporter passed `false` as `PeakDef::peak_json`'s fourth argument;
that parameter is an integer alpha override, so `false` serialized every production Gaussian curve as fully
transparent `rgba(...,0.00)`.  The fitted peaks and continua were present in JSON and tables, but only the
continuum was visible on the chart.

The reporter now inserts `<meta charset="UTF-8">` into the D3 export head before its title/assets.  Production
peak JSON now uses an explicit blue color with alpha override 254, which serializes as 309 opaque
`rgba(47,111,221,1.00)` peak curves in the 33-spectrum smoke gallery; zero matching transparent peak colors
remain.  The regenerated gallery's 37 inline scripts parse successfully and its evaluation objective/counts
are unchanged.

## User-review follow-up: chart defaults and ROI navigation

The gallery now requests log-y unconditionally, including supplied-background reports.  The redundant orange
lower/upper ROI reference lines were removed; ROI extents remain visible through the fitted peak/continuum
shapes and are exact in the diagnostics table.  The derived foreground-minus-scaled-background spectrum is
no longer constructed or plotted by default because the displayed production peaks and continua are fitted
to the original foreground.  The scaled supplied background remains available as context with its explicit
live-time factor.

As a convenience, each chart installs a local click handler on its rendered production `peakFill` and
`peakOutline` paths.  The clicked peak energy is matched to the exact public ROI bounds already stored on the
corresponding diagnostics row; that row is highlighted and scrolled into view.  No D3 library or fit behavior
was changed.  The regenerated smoke gallery contains 33 peak-click handlers, 163 addressable ROI rows, no
background-subtracted curve, no ROI-boundary reference series, and an explicit `setLogY()` call for every
chart.  All 37 inline scripts parse successfully and the objective/counts remain unchanged.
## Manual Detective-X tuning start — 2026-07-20

Starting tree was `2b9abe04` with no tracked diff. All then-untracked artifacts were copied before edits to `/Users/wcjohns/coding/InterSpec_peak_fit_improve/manual_tuning_artifacts/backups/20260720T200322Z-2b9abe04/` (851 MiB; 165 regular-file SHA-256 values verified, plus size and symlink manifests). Focused reporter self-test passed. No tuning had occurred at the time of backup.

Development-only `--spectrum-key-file` now enforces an exact canonical allowlist before fitting, retains canonical corpus order, and rejects duplicate, unknown, and empty selections. `--eval-results-tsv` exports the same cached `ConfigEvaluation` used by the gallery, including objective components, fit state, ROI extrema, continuum types, diagnostics, warnings, and errors. `compare_nuclide_eval.py` strictly joins keys, verifies the per-row objective identity, and emits aggregate, delta, and visual-inspection queue artifacts.

E000 is frozen at `manual_tuning_artifacts/runs/HPGE-DETX-BAL-300-BG-R6OFF-E000/`: Detective-X, Baltimore, 300 s, supplied long background, R6 disabled, no holdout, fixed rel-eff cap 25, eight threads. Diagnostic baseline: 11/11 Success, 0 empty, 0 mechanical failure, total cost -147.686, area cost 0.845927, find reward 7.43726, candidate reward 7.40120, summed miss fraction 0.113311, 6 definite misses, 0 spurious. The usable 11-case sentinel baseline: 11/11 Success, total cost -232.603, area cost 0.961918, find reward 11.2539, candidate reward 11.3083, summed miss fraction 0.0909091, 1 definite miss, 0 spurious. Requested `Co60_Sh` exists as files but is absent from the compact loader's canonical source list; the original 12-key preflight failed before fitting and this is a corpus/metadata defect rather than a silent replacement.

I personally rendered and inspected all 11 diagnostic log-y plots. `Am241_Sh` is bad: an 11-peak shared Quadratic ROI spans 292.10–344.46 keV (52.37 keV, 108 channels, 30.37 FWHM), carries a definite miss, and is visibly far too broad. This is ROI construction/grouping, not an auto-keep-threshold symptom. `Ba133_Unsh` has two definite misses and a large report-only area pull (23.59). `U238_Unsh` has no public ROIs, two definite misses, and an empirical rel-eff fallback warning; classify as low-statistics plus rel-eff/fallback investigation. `U238_Sh` has one definite miss. `Am241_Unsh` also uses fallback and has high SNIP mismatch (7.39), though its FlatStep ROI is visually plausible. Co60 and Cs137 are good ordinary sentinels; Eu152 cases are broadly acceptable dense/multi-ROI examples, with the unshielded low-energy five-peak 12.25-FWHM ROI retained as a watch case. The scalar objective does not adequately expose the Am241_Sh ROI pathology.

Next experiment remains the predeclared E001 auto-keep-significance bracket (2.5 / parent 3.0 / 3.5), on the diagnostic panel first. It tests weak auto-stage U238/Ba133 misses, not the separate Am241_Sh ROI defect. Retention requires a visually defensible recovery with no new sentinel miss/spurious/failure/illegitimate-empty event and no ROI/continuum regression.

### E001–E004 result and structural ROI finding

After extent-aligned scoring was introduced, the default diagnostic was rerun as E000R before comparison.
E001 lower `auto_keep_significance_z=2.5` made a small scalar improvement but recovered no definite
truth and made the bad Am241_Sh representation still more merged; E001 high 3.5 regressed. Both rejected.
E002 `merge_tail_z=3.0` was inert. E003 `initial_nuc_match_cluster_num_sigma=1.0` did not alter Am241_Sh
and only removed a weak U238_Sh fitted peak; rejected. E004 `auto_rel_eff_cluster_num_sigma=1.5` regressed
the diagnostic objective (-147.611 vs E000R -147.782) with no correctness event; rejected.

Source tracing found an existing synthetic-minimum “break too wide cluster” pass, but it runs before later
ROI assembly/refinement. The final public Am241_Sh ROI can therefore remain 30.37 FWHM despite the earlier
cap machinery. This is an implementation invariant gap, not a further threshold-tuning opportunity. A
bounded post-refinement transactional catastrophic-ROI splitter is under focused implementation: it will
only accept a Success result that preserves requested anchors and is no worse in filtered per-channel fit
quality. R4 remains disabled and is not used as production output.

### Final ROI splitter prototype rejected after adversarial review

The first final-width-splitter prototype was compiled and subjected to a focused synthetic helper test, but
was intentionally removed before any corpus evaluation. Independent review found that it mixed fitted and
working calibration frames, selected/split on only fitted peaks while missing unmodeled peak protection,
and compared raw per-channel training chi-square across different channel domains while multiplying
continuum degrees of freedom. Its challenger also included unrelated refinement changes, did not run after
all late ROI additions/overlap resolution, and did not fully protect mixed fixed user ROI geometry. These
are structural correctness defects, not acceptable tradeoffs. The oversized-ROI problem remains open; any
replacement must be a post-final-assembly, solve-calibration proposal that changes only protected eligible
automatic ROIs, rejects unmodeled-peak intersections, and is accepted with a common-domain,
complexity-aware criterion plus end-to-end tests.

### E005 predeclared: manual keep-gate bracket

E005 isolates `manual_keep_significance_z`, whose Detective-X default is 2.0. The lower arm is 1.5 and the
higher arm is 2.5; all other dumped default values remain byte-identical. Unlike E001, this gate operates in
the manual/initial rel-eff stage implicated by the U238 fallback warnings. Hypothesis: 1.5 may retain a
visually defensible weak U238 requested line that is otherwise rejected before auto refinement; 2.5 tests
whether the setting is active and whether marginal manual seeds are driving false structure. Risks are extra
manual candidates, rel-eff contamination, and more/wider ROIs. Acceptance requires recovery of a definite,
visually supported line or another documented manual-stage correction, with no diagnostic/sentinel
miss/spurious/failure/illegitimate-empty regression and no new implausible ROI/continuum behavior. Both arms
will first use the frozen 11-spectrum diagnostic panel and eight threads.

E005 completed through persistent eight-thread terminal sessions; early same-ID directories without a TSV
are preserved as runner-timeout artifacts, not fitter failures. The 2.5 arm produced cost -148.414 versus
-147.782 (area -0.1716, find reward +0.3399, candidate reward +0.1198), but retained all four definite
misses and zero spurious/failure changes. Its meaningful change was Eu152_Unsh (+2 fitted, +1 observable),
not the targeted U238 manual/fallback path. The 1.5 arm regressed to -147.092, again retained all four
misses, and reduced find/candidate reward. Neither arm meets the recovery criterion; do not retain either.

### E006 automatic ROI policy diagnostic regression — classification, 2026-07-21

E006 (`f3ce6c99`) is not corpus-validated. Against E000R on the frozen 11 Detective-X/Baltimore/300 s
diagnostic spectra it regressed from 4 to 11 definite misses (cost -147.782 to -132.075), with 11/11
mechanical Success and no spurious peaks. The net seven-count change is **eight new Am241_Sh misses offset
by recovery of the prior U238_Sh miss**, not seven independent affected spectra. The new Am241_Sh definite
misses are 98.97, 102.98, 146.55, 322.52, 332.35, 376.65, 383.81, and 419.33 keV; 169.56 keV was already
missed in E000R. These lines have strong truth significance (about 8--28 sigma for the low-energy group and
about 6--27 sigma for the 322--419 keV group), so this is not low-statistics ambiguity. The policy does
reduce the previous 30.37-FWHM / 11-peak Am ROI to smaller ROIs, but the quality loss is unacceptable.

Source trace classifies the dominant cause as a **policy-integration/model-evidence propagation defect**,
not a default-weight choice. In `cluster_gammas_to_rois`, the policy-separated path chooses a boundary using
only each adjacent group's largest-amplitude gamma as its anchor. Its `left_valid` / `right_valid` checks
then verify only those dominant anchors. A boundary can therefore retain those two anchors while clipping
other requested gammas belonging to either multi-line group; the later under-width path can discard the
clipped automatic cluster. The existing non-policy fold-back guard is explicitly disabled when
`use_automatic_roi_policy` is true, so it cannot preserve the previously admitted gamma. A correct policy
split must preserve *every* requested gamma in both groups (or decline that split/repartition before any
geometry is committed); it must never use dropping an admitted group as a way to satisfy separation.

The U238_Sh recovery is a separate favorable result but does not compensate for the Am correctness failure.
Do not tune policy weights until strict full-group preservation is added, focused-tested, and E006 is rerun.

### Planned final-ROI guard redesign (not implemented)

**Purpose.** Enforce a final automatic-ROI width rail after all ordinary refinement, escape/511 additions,
and overlap resolution have assembled the geometry actually passed to the solver. It is not R4 and must not
become a general boundary optimizer.

**Eligibility.** Consider only a final automatic ROI whose width exceeds
`auto_rel_eff_sol_max_fwhm` in the current solve calibration. Exclude every user-supplied fixed ROI, every
ROI that intersects protected mixed/user geometry, and an R2-rescue iteration. Derive the protected ranges
from all user geometry, not only the current `existing_roi_ranges` subset.

**Proposal.** Work in the same energy calibration as the foreground supplied to the final solve. Start from
the incumbent final ROIs and replace only one eligible oversized parent at a time. Candidate boundaries must
be clean windows between adjacent requested/fitted anchors, must not intersect an unfit/auto-peak core or
tail, and must have an observed valley compatible with the shared global continuum. Children retain complete
coverage of the original parent; no channel may be silently dropped. Use the fewest boundaries that bring
each child under the cap, otherwise leave the parent unchanged.

**Transaction.** Solve the incumbent and challenger on one identical, immutable channel domain. Compare a
complexity-aware fit criterion (AICc/BIC or an explicitly justified likelihood penalty) that includes all
new per-child continuum parameters. Reject if the challenger fails, loses any requested source anchor,
creates a new bad overlap, worsens the criterion, or its *actual final solved* ROI geometry still violates
the claimed split. Do not allow unrelated refinement, continuum, or peak-set changes to supply the
acceptance improvement.

**Tests before corpus use.** Add end-to-end fixtures for: calibration refinement; an unfit peak in a would-be
gap; fixed/mixed user ROI protection; escape/511 re-merge after the guard; strict requested-anchor
preservation; common-domain/complexity-aware acceptance; equal-quality but width-correct geometry; and the
actual one-channel gap in solve calibration. Then run both builds, reporter tests, the full source-fit test
suite, diagnostic and sentinel A/B galleries.

## Central automatic ROI merge/partition policy — structural; not corpus-validated

### Architecture

This change replaces the proposed final-width splitter with one pairwise automatic-boundary policy.  The
policy consumes adjacent groups with their modeled energies/areas, the solve's current foreground and
calibration-matched global SNIP estimate, current-frame unfit peak cores, the fitted FWHM function, and the
existing merge/width/continuum penalties.  It emits a typed `AutomaticRoiDecisionDiagnostic` through
`PeakFitResult::automatic_roi_diagnostics`, including separation/FWHM, observed/SNIP valley counts, modeled
tail and unexplained-excess significance, SNIP mismatch, child sideband channels, width pressure, identical-
coverage one/two-continuum AICc scores, calibration channel count, boundary channel, decision, and reason.

Candidate boundaries are channel-aligned and require at least eight usable sideband channels per child.
The comparison fits the existing Linear, Quadratic, FlatStep, and LinearStep families to the shared SNIP
continuum with raw foreground Poisson variances and rejects negative/non-finite predictions.  One and two
continua are scored on the same union of channels with the configured AICc penalty; the one-ROI score also
receives the soft quadratic width pressure.  A clean boundary must also pass the existing tail and gap
thresholds and unexplained-excess test.  An unmodeled +/-1-FWHM core is an exclusion, never a legal split or
merge-through.  Wide joins are allowed only as the explicit `MergeInseparableWide` outcome.

FitPeaksForNuclides-owned automatic ranges are finalized as `Fixed` before `RelActCalcAuto::solve`, closing
its generic `CanBeBrokenUp` recombination bypass.  Mixed fixed user ranges are tracked separately, restored
by original channel coordinate after each calibration update, and not resized or merged.  Unfit peaks and
the SNIP estimate are likewise re-materialized/re-stamped in the current working calibration.  The policy
scores continuum plausibility only; it does not replace the ROI's existing selected continuum family.

### Routed sites and explicit bypasses

In policy mode the decision is used by `cluster_gammas_to_rois` (manual, fallback, and ordinary refinement),
`merge_rois` (top-level source/NORM combination and the in-solve NORM combination), initial/no-ecal/
desperation escape-overlap finalization, refinement edge/found reconciliation, refinement desperation
finalization, and refined-solve escape-overlap finalization.  Tight found-peak and escape insertion remain
simple constructors, but every overlap they create is reconciled by the policy before solve.  Initial
clustering now returns its modeled energy/area vectors to late reconciliation; refinement carries those
vectors plus current fitted predictions, while floating escape peaks use matching observed peak area when
available.  Missing model evidence rejects the late automatic addition instead of inventing a midpoint core.
The old unconditional refinement edge-fold loop is gone.  Channel-disjoint ranges are made by the policy
resolver; `ensure_min_channel_gap` remains only on legacy branches.

The explicit R6-enabled mode is fenced to the pre-policy path end-to-end: manual/fallback/no-match
construction, source-only incumbent geometry, tight-seed range types, overlap resolution, and RelActAuto
recombination retain their prior inputs and behavior.  Only an actually accepted R6 selected-line coverage
union or retained nuisance neighborhood emits `R6LegacyBypass`.  R2's private reclassification/rebuild
clustering and non-overlapping transactional admission remain legacy and unchanged.  511 keV remains
floating-peak-only; R4 remains shadow-only and disabled.

### Parameters and defaults

No gene, serialized field, command-line option, or default numeric value was added or changed.  Existing
`merge_tail_z`, `merge_clean_gap_fwhm`, `cont_order_aicc_penalty`, `auto_roi_core_num_fwhm`, and
`auto_rel_eff_sol_max_fwhm` now influence more automatic joins.  The last value is a soft pressure onset, not
a hard cap.  Detector-class manual re-evaluation is therefore required even though defaults are byte-stable.

### Verification status

Both required builds completed successfully:

- `cmake --build target/peak_fit_improve/build_vscode`
- `cmake --build target/testing/build_test --target test_fitPeaksForSources`

Focused `StatisticalDetailHelpers` completed 9/9 with no errors, including deterministic clean-valley,
overlap, mean-outside/core-crossing unmodeled protection, unavailable-AICc, invalid-resolution,
shifted-calibration, protected-geometry, and sequential 11-group width-pressure cases.  The reporter
self-test passed 1/1.  The final rebuilt full `test_fitPeaksForSources` execution passed all 34 cases with
exit code 0.  The R6 default-sequence isolation sentinel passed 248/248 assertions after the final legacy
metadata fence.  `git diff --check` passed.  No GA or corpus evaluation was run.

### Diff summary and risks

The structural diff is limited to `InterSpec/FitPeaksForNuclides.h`, `src/FitPeaksForNuclides.cpp`, focused
tests in `target/testing/test_fitPeaksForSources.cpp`, and this review record.  Dirty tuning-harness files and
all other tracked/untracked artifacts were preserved.  No commit was made.

Remaining risks are that joined-group lineage resets across separate late resolver invocations; a modeled
line with no fitted/observed area has a conservative zero-tail estimate; unmodeled-core exclusion currently
rejects the whole anchor gap rather than searching around multiple occupied candidates; and no deterministic
fixture directly forces each of found, restored-edge, and escape insertion (the full-fit invariant instead
checks that policy-mode solve ranges are channel-disjoint and `Fixed`).  The runtime-only legacy splitter is
retained solely inside the explicit R6 compatibility fence.  These are structural follow-ups for corpus
evidence, not authority to weaken the clean-gap or protected-geometry invariants.

### Exact first primary-session corpus experiment (only after blockers are green)

Run the frozen 11-spectrum Detective-X/Baltimore/300-s diagnostic panel, supplied background, R6 disabled,
and E000R genes:

```sh
mkdir -p /Users/wcjohns/coding/InterSpec_peak_fit_improve/manual_tuning_artifacts/runs/HPGE-DETX-BAL-300-BG-R6OFF-E006-ROI-POLICY/diagnostic
cd /Users/wcjohns/coding/InterSpec_peak_fit_improve/manual_tuning_artifacts/runs/HPGE-DETX-BAL-300-BG-R6OFF-E006-ROI-POLICY/diagnostic

/Users/wcjohns/coding/InterSpec_peaks_for_source_opt/target/peak_fit_improve/build_vscode/PeakFitImprove \
  --action NuclideConfigEval \
  --data-base-dir /Users/wcjohns/coding/InterSpec_peak_fit_improve/peak_fit_accuracy_inject_compact \
  --static-data-dir /Users/wcjohns/coding/InterSpec_peaks_for_source_opt/data \
  --det-type HPGe \
  --detector Detective-X \
  --city Baltimore \
  --live-time 300_seconds \
  --spectrum-key-file /Users/wcjohns/coding/InterSpec_peaks_for_source_opt/target/peak_fit_improve/manual_tuning_panels/diagnostic.keys \
  --config-genes /Users/wcjohns/coding/InterSpec_peak_fit_improve/manual_tuning_artifacts/runs/HPGE-DETX-BAL-300-BG-R6OFF-E000R-EXTENT/diagnostic/HPGe_default_genes.tsv \
  --holdout-frac 0 \
  --holdout-role all \
  --number-threads 8 \
  --rel-eff-chi2-cap-mode fixed \
  --rel-eff-chi2-cap-value 25 \
  --disable-auto-interferer-fit \
  --eval-results-tsv results.tsv \
  --eval-html gallery.html

python3 /Users/wcjohns/coding/InterSpec_peaks_for_source_opt/target/peak_fit_improve/compare_nuclide_eval.py \
  /Users/wcjohns/coding/InterSpec_peak_fit_improve/manual_tuning_artifacts/runs/HPGE-DETX-BAL-300-BG-R6OFF-E000R-EXTENT/diagnostic/results.tsv \
  results.tsv \
  comparison \
  --target-key-file /Users/wcjohns/coding/InterSpec_peaks_for_source_opt/target/peak_fit_improve/manual_tuning_panels/diagnostic.keys
```

Acceptance is 11/11 successes, no new empty/failure/spurious event, no increase beyond E000R's four definite
misses, manual inspection of all plots, and removal of the silent Am241_Sh 292--344-keV/30.37-FWHM ROI.  Any
remaining wide ROI must carry a convincing `MergeInseparableWide` diagnostic and be visually defensible
before proceeding to the sentinel panel.

---

# 2026-07-21 — Atom-safe automatic ROI split/combine layer (replaces the E006 policy; NOT corpus-validated)

Supersedes the E006 regression classified above. All changes in the working tree, NOT committed. The
uncommitted 38-line single-site preservation patch was reverted first (saved to scratch) and superseded.

## Architecture

Every automatic (policy-mode) ROI split/combine now operates on stable **atoms** — one per admitted
modeled gamma line (or found/floating feature) — carried WITH the geometry instead of reconstructed from
a flat `(energy,area)` list by geometric containment. New `detail::` layer (header
`InterSpec/FitPeaksForNuclides.h`, impl `src/FitPeaksForNuclides.cpp`):

- `RoiAtom{id,energy,area,kind,source,rel_eff_curve_index,admission}` + `next_roi_atom_id()` (atomic).
- `AutomaticRoiComponent` — channel-aligned bounds + its exactly-once atoms + pass-through
  continuum/range/protected metadata.
- `partition_automatic_roi_pair(...)` — materializes one adjacent pair into channel-aligned children
  (spatial atom assignment) or a merge, using the UNCHANGED `evaluate_automatic_roi_boundary` oracle for
  the decision. Owns: core-safe boundary-channel search (a boundary may never cut an atom's
  ±peak_core_num_fwhm·FWHM core), min-width outward widening, channel-aligned unmodeled-exclusion carves
  (declined when they would clip an admitted core), and protected-edge pinning. No core-safe partition ⇒
  MERGE (or, for a protected pin, retain) — never a dropped side.
- `reconcile_automatic_components(...)` — the single left-fold driver over energy-sorted (possibly
  overlapping) components; folds each adjacent pair, re-examines an enlarged component after a merge;
  validates the whole-stage transaction, works all-or-nothing on a copy (valid==false ⇒ caller retains
  incumbent).
- `validate_automatic_roi_transaction(...)` — always-run: atom-ID multiset preserved exactly-once
  (before == after ∪ orphans), sorted channel-disjoint, atom energy + clamped-core containment (core
  containment is dev-hard / release-soft, since materialized children satisfy it by construction and only
  a pass-through input narrower than its own core can trip it — energy containment, the loss-safety
  property, always holds), bit-identical protected bounds/metadata.
- `assign_atoms_to_disjoint_rois(...)` — exact-once snapshot assignment for a disjoint-ROI bridge.

Additive-only diagnostics: `AutomaticRoiDecision::InfeasiblePartition` (appended);
`AutomaticRoiDecisionDiagnostic{atoms_reassigned, partition_infeasible}`. NO gene / serialized field /
CLI option / default numeric value changed. New `GammaClusteringSettings::cluster_admission_stage`
(diagnostics-only enum, not serialized/tuned).

## Routed policy call sites (all three)

1. `cluster_gammas_to_rois` — atoms minted at the keep-gate accept (and on `MarginalRejectedCluster` for
   R2), carried on `ClusteredGammaInfo`. Policy-mode merge loop calls `partition_automatic_roi_pair` per
   adjacent pair (then `continue`), bypassing the prior split/valley/carve code and its
   pop_back/continue drop paths. The under-width emit path folds-with-atoms into an abutting previous ROI
   or emits the narrow ROI — never `continue`-drops admitted atoms in policy mode.
2. `merge_rois` — a policy branch reconciles all ROIs at once (atoms minted from each ROI's modeled
   lines), replacing the per-pair split/merge with its pop_back+reprocess / skip-current drop paths.
3. `resolve_automatic_overlapping_rois` — a policy path assigns each modeled candidate and floating
   feature EXACTLY ONCE to its containing ROI (nearest-midpoint on ties, replacing the inclusive
   containment scan that double-claimed atoms across overlapping ROIs), then reconciles; the
   insufficient-evidence / invalid-child continue drop paths are gone. Honors the reconcile `valid`
   verdict: on failure it retains incumbent geometry via the legacy channel-disjoint resolver.

## Preserved / explicit remaining bypasses

- `use_automatic_roi_policy == false` (R6-enabled fits): byte-identical legacy split/merge and legacy
  `resolve_overlapping_rois`+`ensure_min_channel_gap` in all three functions; also the no-calibration
  fallback for policy mode.
- R2 rescue admission rule + its private `use_automatic_roi_policy=false` re-clustering: unchanged.
- Escape/511/found-peak-seed insertion + admission: unchanged; they participate as floating/evidence
  atoms during resolve exactly as their flat-list evidence did before.
- R4 shadow optimizer: still disabled/diagnostics-only.

## Exact-once preservation argument

Within each policy operation the atom-ID multiset is conserved: `partition_automatic_roi_pair`
concatenates all atoms on merge, and on split assigns every union atom to exactly one child by energy
across a gap that (by the core-safe test) contains no atom; orphans arise only from a protected-boundary
straddle and are reported, never silently discarded. The reconcile driver validates
`before == after ∪ orphans` (each id once) and rejects any non-disjoint result; the zero-atom collision
branch only ever discards an empty, non-protected side (the non-empty side is always kept), and defers a
protected-empty collision to the protected pin. Dev builds assert on any invalid transaction; release
degrades to the merge fallback or retains incumbent — never a crash.

## Adversarial review (subagent) — 2 findings, both fixed

1. HIGH — the reconcile zero-atom branch used an OR guard that could `continue`-drop the NON-empty side
   when the wider ROI was fully starved of atoms by a nearer overlapping ROI (a real exposure of the new
   exact-once nearest-midpoint assignment); callers also ignored the `valid` verdict. FIXED: the branch
   now discards only the empty non-protected side and keeps the atom-bearing side; `resolve` honors
   `rr.valid` (legacy-resolver fallback). Regression test `test_reconcile_starved_wide_roi_keeps_atom`.
2. LOW — `atomlayer_partition_protected` at a protected ROI whose lower edge is spectrum channel 0 could
   emit a channel-0 overlap. FIXED: that degenerate case now dissolves the (zero-room) automatic span.

## Tests

Ten new deterministic `AutomaticRoiPartition`/reconcile unit tests (in `StatisticalDetailHelpers`), each
asserting no atom loss/duplication via a shared `check_exact_once`: spatial reassignment of a
non-dominant multi-line member (the Am241 shape), channel-rounding boundary, no-core-safe-channel→merge,
under-width child at the spectrum edge, exclusion-band-clips-core declined, protected-geometry pinning +
straddle-orphan, overlapping-ROI exact-once ownership, the starved-wide-ROI regression, evidence-only
survival, and a seeded 300-iteration randomized exact-once sweep. `verify_fit_result` additionally
asserts every `partition_infeasible` diagnostic carries a reason (plus the existing Fixed +
channel-disjoint policy-ROI checks).

## Pragmatic scope note (residual limitation)

`merge_rois` and `resolve` mint fresh atoms per call from their existing modeled inputs rather than
threading atom IDENTITY end-to-end through `fit_peaks_for_nuclide_relactauto`'s signature. Each policy
operation is therefore independently exact-once (satisfying the required no-loss/no-dup/no-misassign
invariant), but cross-stage provenance for equal-energy lines from different sources is not carried
between stages. The top-level source+NORM `InitialRoi` build still associates the flat candidate list by
containment — exact-once over the disjoint cluster/merge output, but could double-count a candidate lying
in an overlap between a source ROI and a NORM ROI before `merge_rois` folds them (a pre-existing edge
case, not the E006 regression). Full identity threading is a deferred follow-up.

## Verification status

- `cmake --build target/peak_fit_improve/build_vscode` — clean.
- `cmake --build target/testing/build_test --target test_fitPeaksForSources` — clean.
- `ctest -R '^PeakFitImproveReporterSelfTest$'` — 1/1 passed.
- Full `test_fitPeaksForSources` (PERFORM_DEVELOPER_CHECKS=ON) — 44/44 cases, "No errors detected".
- `git diff --check` — clean. No GA or corpus evaluation run.

Diff summary: `InterSpec/FitPeaksForNuclides.h` (+170), `src/FitPeaksForNuclides.cpp` (~+1400/-203),
`target/testing/test_fitPeaksForSources.cpp` (+520). Tuning-harness + untracked artifacts preserved; no
commit.

**Per-detector-class manual re-evaluation + user visual review are required before these policy-mode
geometries are trusted** — the partition materializes boundaries differently (core-safe channel vs
valley split) and now merges-rather-than-drops where the old policy dropped. GA remains deferred.

## Risks

- Policy-mode ROI geometry changes intentionally vs `f3ce6c99`; only the user-owned frozen E006 rerun
  confirms Am241_Sh recovery.
- `merge_rois` min-width expressed at the child midpoint vs the old per-side `width>=FWHM` (policy-only,
  semantically equivalent, not bit-equal).
- Cross-stage atom identity not threaded (see scope note).

## Exact frozen E006 command the primary session should run next

The E006 command block recorded above (Detective-X / Baltimore / 300 s diagnostic panel, supplied
background, R6 disabled, E000R genes) is unchanged; rerun it and diff against E000R with
`compare_nuclide_eval.py`. Acceptance: 11/11 Success, ≤4 definite misses (recover the eight Am241_Sh
98.97–419.33 keV misses), no new spurious/failure/empty, and the silent 30.37-FWHM Am241_Sh ROI gone or
carrying a visually defensible `MergeInseparableWide` diagnostic.

### E006 rerun with the atom-safe layer (2026-07-21) — NO recovery; regression is a SPLIT-DECISION problem, not atom-safety

Ran the exact frozen E006 command above (atom-safe binary, R6 off, E000R genes, 11-spectrum
Detective-X/Baltimore/300 s panel). Result is essentially IDENTICAL to the buggy `f3ce6c99` E006:

| metric | E000R (parent) | buggy E006 (`f3ce6c99`) | atom-safe |
|---|---:|---:|---:|
| scalar_cost | -147.782 | -132.075 | **-133.373** |
| definite_miss_count | 4 | 11 | **11** |
| success / spurious | 11 / 0 | 11 / 0 | 11 / 0 |

Per-spectrum: the net is unchanged — **Am241_Sh Δmiss = +8**, U238_Sh Δmiss = -1, all others 0. On
Am241_Sh the atom-safe layer DID shrink the pathological ROI (max ROI 30.4→15.6 FWHM, 52→30 keV,
11→6 peaks/ROI) but that shrink is exactly what loses the lines: observable peaks 23→14, fitted peaks
52→26, definite misses 1→9.

**Diagnosis — the E006 root-cause classification (atom drops in `cluster_gammas_to_rois`) was
INCOMPLETE.** The atom-safe layer provably keeps every admitted line in a channel-disjoint ROI (the
invariant holds; full suite green with dev-checks ON), yet the Am241 misses persist. The lines are NOT
dropped from the ROI set — they are pruned by the FIT/observable filter once the wide multi-line Am241
ROI is split into narrower ones. The pathological 30-FWHM Am241 ROI, while visually ugly, was
FUNCTIONALLY NECESSARY: co-fitting its ~11 lines on one shared continuum + rel-eff curve is what made
the weak 98–419 keV lines observable. The policy's decision to split it (the oracle sees clean valleys
between the well-separated HPGe lines and returns KeepSeparate) is the real regression driver, and the
atom-safe layer faithfully executes that split without dropping atoms — same misses.

**Consequence / next direction.** Atom-safety is necessary but NOT sufficient to recover E006. The fix
to pursue is the split DECISION for weak multi-line sources: the policy must recognize that splitting a
group whose weak members are only findable via the shared rel-eff/continuum is harmful, and keep it
merged (or the observable-peak filter must credit rel-eff-anchored weak lines the way the wide ROI did).
This is model/decision work on `evaluate_automatic_roi_boundary` (or the post-fit observability test),
not further atom-safety plumbing. The one genuine improvement here is U238_Sh (-1 miss, +1 observable).

Artifacts (temp): `/private/tmp/E006-atom-safe-rerun/diagnostic/` — `results.tsv`, `gallery.html`
(4.96 MB SpectrumChartD3 gallery), `comparison/` (aggregate + per-spectrum delta vs E000R). The
`resolve_overlapping_rois` overlap warning at 1457–1464 keV is pre-existing (present in the buggy E006
`developer_errors.log`), not introduced by this change.

### Am241_Sh miss mechanism traced (2026-07-21) — CORRECTS the direction above: policy OVER-MERGES, collapsing the activity solve

Traced the three user-flagged strong lines (322.5, 376.6, 383.8 keV; truth 16–20 sigma) with
`NuclideFitDebug --det HPGe --src Am241 --disable-interferer --trace` (added a `--disable-interferer`
flag to that untracked tool to reproduce the eval's policy mode). The debug run reproduces the eval
exactly (observable=14, fit=26).

**They are NOT lost to any atom-drop path** — all pass the initial keep-gate at z=11–14 and are in ROIs
at initial clustering. The exact exclusion point is the **refinement re-clustering keep-gate**
(`cluster_gammas_to_rois`, the `est_counts > sm_keep_gate_min_est_counts`=15 floor): after the first
solve the fitted Am241 activity **collapses to 16.6M Bq** (correct value ~130M), so the refinement
predicts sub-floor counts (376.6: est_counts 363→3.5; 322.5→1.0; 383.8→0.8; 332.4→1.3) and rejects them;
the final solve then has no peak for them at all.

**Root cause = policy OVER-MERGE, not a split.** Initial ROIs, same binary, same spectrum:
- Policy merges 293–348 keV into ONE 30.1-FWHM Quadratic ROI (17 gammas) and 361–391 into one 15.5-FWHM
  ROI. That wide flexible ROI collapses the RelActCalcAuto activity solve.
- Legacy (R6 path) splits the same region into FIVE narrow ROIs (6–7 FWHM) + three at 361–391, keeps the
  activity STABLE at 131→126→127M, and finds all seven lines with correct areas (observable=23).

The eval's "policy shrank the ROI to 15.6 FWHM" is the DAMAGED post-collapse remnant (refinement dropped
most lines, shrinking the ROI), NOT a healthy split. This corrects the earlier note that read the wide
ROI as "functionally necessary" — it is the opposite: the wide ROI is what breaks the fit.

**Mechanism of the over-merge.** `f3ce6c99` disabled the `break-too-wide` splitter in policy mode
(`final_clusters = merged_clusters`), relying on the oracle's `max_width_fwhm` as a SOFT width pressure.
For the dense Am241 forest the clean-gap test fails (weak lines contaminate the gaps) so the oracle
MERGES, and with no hard width cap the ROI grows to 30 FWHM. **Fix direction: a hard width cap /
`break-too-wide` equivalent in policy mode (the "final-ROI width guard" already sketched above), so the
forest is split into solve-stable narrow ROIs.** The atom-safe layer is the correct prerequisite (it
guarantees the split won't drop lines); the width cap is the missing decision piece.

## Shielded-source automatic recovery experiments and retained change — 2026-07-21

### Corrected boundary semantics: absence of a significant peak bridge

The phrase “credible clean valley” was too morphological.  A noisy or low-statistics spectrum need not
show a literal local minimum when the peaks on both sides are insignificant.  The boundary evidence is
now documented and diagnosed as **no significant peak bridge**: predicted Gaussian-tail content from
both modeled groups must be insignificant over the proposed boundary window, sidebands must support the
shared-continuum estimate, and there must be no statistically significant positive unexplained
peak-like excess over it.  This is positive evidence that significant peak counts do not connect the
ROIs; it is not a demand for a low raw-count trough.  The existing statistics and thresholds are
unchanged; only the interpretation/name was corrected.

This also corrects the causal framing of the 52.45-keV / 30.32-FWHM Am241_Sh union.  Shielding does not
create that low-energy union.  The seed/merge/boundary mechanism fails to repair it.  Shielding makes the
bad geometry much more damaging because the empirical rel-eff/activity solve is ill-conditioned and the
over-wide flexible continuum can swallow the high-energy anchors, causing the later 15-count keep-gate
collapse.

### Frozen evidence

- Pre-change working-tree patch: `/private/tmp/shielded-source-recovery-baseline/prechange-working-tree.diff`.
- Frozen comparable Detective-X/Baltimore/300-s corpus (205 spectra):
  `/private/tmp/shielded-source-recovery-baseline/full_baltimore_300/results.tsv`: 80 definite misses,
  7 spurious, 0 mechanical failures, 5 legitimate empties, scalar cost -2860.70220327.
- Frozen atom-safe 11-spectrum panel:
  `/private/tmp/E006-atom-safe-rerun/diagnostic/results.tsv`: 11 definite misses, 0 spurious.
- Legacy Am241 trace stayed at 131.1 -> 126.2 -> 127.5 MBq with 23 observable / 51 fitted peaks.
  The final legacy trace has the same behavior.  The observable N42s are identical after ignoring only
  `n42DocDateTime`; the private fit N42 writer changes continuum/peak serialization order, but a
  structure-aware comparison found identical multisets of all 20 continua and all 305 peaks.  Thus the
  requested byte comparison is semantically exact; literal bytes are precluded only by timestamp and
  nondeterministic serialization order, not a fit difference.

### Experiment ladder

1. **External-attenuation physical manual candidate — rejected.**  The fixed-Z shielding candidate
   reused the configured shielding material, areal-density bounds, Hoerl eligibility, and AICc parameter
   accounting.  On Am241_Sh it drove the activity to about 4.17 GBq (rather than the 126–131 MBq legacy
   reference) with attenuation/Hoerl terms pinned at bounds.  It did not establish a physically stable
   root fix, so all production shielding-candidate code was removed.  No fit-Z parameter, GA gene, or
   source-specific switch was added.
2. **Whole-union scored geometry challenger — rejected.**  Width pressure was evaluated before the
   overlapping-anchor short circuit and core-safe partitions were scored with the existing global-
   continuum AICc plus soft width term.  It helped the diagnostic panel, but the 205-spectrum corpus
   introduced broader regressions.  The experiment, its tests, and its materializer/scorer changes were
   reverted; the authoritative atom-safe geometry is otherwise unchanged.
3. **Internal over-wide single-component splitter — rejected.**  A second geometry trial also recovered
   the exemplar but did not pass the broad gates.  It was fully removed; the legacy hard splitter and a
   new width constant were not introduced.
4. **Broad all-significant-anchor refinement protection — rejected.**  Keeping all significant anchors
   through refinement produced a new Pu238_Sh spurious peak and a U235 regression.  This demonstrated
   why blanket found-peak protection is unsafe.  The experiment was removed rather than stacked.
5. **Source-clean transactional seed challenger — retained.**  This is the smallest candidate that
   passes the panel and improves the comparable full corpus.  The manual stage carries a separate,
   contaminant-pruned set of data-significant, FWHM-distinct requested-source anchors and source-clean
   alternate ROIs.  Normal initial ROIs still contain the original found/matched seeds.

### Retained automatic-only transaction

After the first empirical automatic solve, recovery is attempted only when at least two independent
clean requested-source anchors that were observed significantly would fall below the existing strict
`predicted counts > 15` refinement gate.  The clean alternate initial ROIs are then solved as a
transactional challenger.  It replaces the incumbent only when all of the following hold:

- the solve succeeds and preserves strictly more predicted clean anchors above the existing keep gate;
- it preserves at least as many FWHM-distinct significant fitted requested-source peaks, preventing one
  fitted peak from satisfying multiple close anchors; and
- it improves data-only AICc on identical immutable +/-0.5-FWHM windows around the same clean anchors.
  The score uses Poisson deviance and obtains the effective parameter count from
  `solution_data_row_count - m_dof_data`, with the existing manual rel-eff AICc penalty.  A finite
  challenger may replace an incumbent for which the common-channel score is unavailable; two
  unavailable scores reject.

Rejected challenger diagnostics are rolled back to a vector checkpoint.  Once accepted, both recovered
predicted-anchor and fitted-source-evidence floors remain in force for the no-energy-cal/physical retries
and ordinary refinement, so a later empirical solve cannot immediately undo the recovery.  The code is
inside `use_automatic_roi_policy == true`; R6 legacy, R2 admission, protected/user/mixed geometry, NORM
geometry, and atom ownership are unchanged.  There is no production API/serialization change and no new
tunable threshold.

Focused coverage exercises collapse triggering (two losses required), non-collapse/non-shielded
rejection, failed-solve rejection, strict anchor improvement, fitted-source-evidence preservation,
finite/unavailable AICc behavior, and the data-only AICc parameter calculation.

### Mechanistic and corpus evidence

Final Am241_Sh trace (`/private/tmp/aicc_am.{out,err}`): 17 clean anchors; predicted anchors 6 -> 17;
FWHM-distinct fitted source evidence 39 -> 145; common-anchor data-only AICc 7833.7 -> 567.6; challenger
accepted.  The accepted activity is 108.486 MBq, within a factor of two of the 126–131 MBq legacy
reference instead of the original 16.6-MBq collapse.  The eight missed lines are observable at
98.915, 102.943, 146.381, 322.577, 332.290, 376.673, 383.931, and 419.101 keV; in particular the
322/376/383-keV anchors remain above re-clustering and survive observable refit.  Final result: 25
observable / 59 fitted peaks.  A U238_Sh guard trace had two anchors but improved only 0 -> 1 predicted
anchor while reducing fitted evidence 11 -> 6 and having unavailable common scores, so it was correctly
rejected.

Final sentinel traces on the identity-guard build also completed successfully without invoking the
source-clean challenger: Am241_Unsh (1 observable / 1 fitted), Eu152_Sh (17 / 27), Eu152_Unsh
(29 / 48), and Cs137_Unsh (1 / 1).  Their traces are `/private/tmp/identity_<source>.{out,err}`.

Final 11-spectrum panel (`/private/tmp/shielded-source-recovery-retained3/diagnostic`): 11/11 Success,
3 definite misses, 0 spurious, 0 failure, scalar cost -150.520649585.  Versus atom-safe E006 this is
-8 misses, all from Am241_Sh, with no new miss/spurious on the other ten.  Versus E000R it is one fewer
miss (3 vs 4), no spurious/failure, and a 2.738 lower scalar cost.  The U238_Sh improvement is retained.
Comparisons are in `comparison_atomsafe_final/` and `comparison_e000r_final/`; the complete gallery is
`gallery.html`.

Final comparable 205-spectrum corpus
(`/private/tmp/shielded-source-recovery-retained3/full_baltimore_300`): definite misses 80 -> 71,
spurious 7 -> 7, mechanical failures 0 -> 0, legitimate empties 5 -> 5, and scalar cost
-2860.70220327 -> -2878.16818763 (improvement 17.46598435).  Miss changes are Am241_Sh -8,
U235_Unsh_9000 -2, and U235_Unsh_8000 +1; there are no spurious-count changes.  This therefore improves
beyond the single shielded exemplar by eliminating the same anchor-collapse class in U235_Unsh_9000.
The comparison and full gallery are `comparison/aggregate_comparison.md` and `gallery.html`.

The compact external corpus now enumerates 1,819 usable Detective-X spectra across additional cities
and live times; the frozen apples-to-apples “full” gate available for this task is the established 205
Detective-X/Baltimore/300-s set.  An accidental 1,819-spectrum precompute was stopped before evaluation
and produced no comparison result.  The stable 205-spectrum gallery, not a partial expanded run, is the
evidence cited here.

### Adversarial review and residual risk

The requested independent adversarial review found five substantive issues, all addressed before the
final panel/corpus reruns:

1. later retries/refinement could undo an accepted recovery — fixed with persistent predicted/fitted
   evidence floors;
2. the first AICc comparison used each solution's own channels — fixed with identical anchor-local
   channels and parameter counts from data rows/`m_dof_data`;
3. one fitted peak could appear to support multiple close anchors — fixed by counting FWHM-distinct
   requested fitted evidence exactly once; and
4. rejected challenger diagnostics leaked — fixed with checkpoint/resize rollback; and
5. downstream floors initially preserved only cardinality, allowing one recovered line to be exchanged
   for another — fixed by retaining and checking every above-gate anchor and every FWHM-distinct fitted
   requested gamma by identity.  The final identity-guard panel and corpus TSVs are byte-identical to
   the preceding accepted run.

The final reviewer re-check found no remaining HIGH or MEDIUM correctness/isolation issue.  A LOW
theoretical caveat remains for fitted lines without an assigned theoretical gamma: the identity guard
uses all-of/any-of rather than consumed-candidate matching, so an unusually large change in the fitted
resolution could in principle let one unassigned line match two required identities.  Assigned gamma
energies and the FWHM-distinct construction prevent this for the normal requested-nuclide path.  The
focused selection tests do not inject an end-to-end throwing challenger or an identity-swap refinement;
the transactional rollback and identity guards are covered by review plus the byte-identical panel/full
reruns rather than a dedicated forced-failure fixture.

Residuals requiring visual/user review: U235_Unsh_8000 gains one definite miss while U235_Unsh_9000
loses two; the largest scalar regressions are U235_Unsh_2000 +2.794, U235_Unsh_8000 +2.013, and W187_Sh
+1.356.  Some very wide Am241 ROIs remain because geometry experiments failed the full gate; the
retained fix prevents their flexible continua from erasing source evidence but does not claim to solve
all ROI aesthetics.  The private fit N42 serializer remains order-nondeterministic.  The implementation
retains the existing 15-count gate and AICc penalty, so detector classes beyond this HPGe corpus still
need their normal per-class manual evaluation.

Final verification on the retained identity-guard build:

- both focused CMake targets rebuilt successfully;
- `test_fitPeaksForSources`: 45/45 cases (the prior 44 plus the new source-clean/AICc case),
  `*** No errors detected` (`/private/tmp/shielded_source_final_tests2.out`);
- `PeakFitImproveReporterSelfTest`: 1/1 passed;
- the final panel and 205-spectrum result TSVs are byte-identical to the pre-identity-guard accepted
  run, while the legacy observable N42 is identical after timestamp normalization and its private fit
  contains the same 20-continuum / 305-peak structural multiset; and
- `git diff --check` is clean.  Nothing was staged or committed.

## E007 — source-evidence-gated component ROI recovery (2026-07-21)

### Why the retained root fix was not enough

The retained source-clean challenger fixed the rel-eff/activity collapse, but its Am241_Sh result still
fit the unsupported 300.46/312-keV predictions in one 293--348-keV ROI with the real 322/332/335-keV
structure.  That remaining defect is geometric.  Shielding does **not** create the 52.45-keV-wide /
30.32-FWHM union; the seed/merge/boundary mechanism fails to repair it.  Shielding only makes the bad
union more damaging by destabilizing the rel-eff/activity solve and allowing its continuum to swallow
real high-energy anchors.

The boundary rule is consequently expressed as **no significant peak bridge**, not a literal clean
valley.  A boundary is supported when modeled tails and positive unexplained peak-like counts are both
insignificant on the boundary channels, with adequate continuum sidebands.  Low-statistics data need not
show a morphological trough to establish that significant peak counts do not connect two ROIs.

### Retained evidence and geometry transactions

The following work remains gated behind an already accepted multi-anchor source-collapse recovery.  It
does not run in ordinary initial policy fits or the R6/legacy path.

1. Each provisional source group in an over-wide overlap component is compared on identical measured
   foreground channels under H0 (production continuum only), Hs (the exact fixed-ratio requested-source
   line mixture with one nonnegative common scale), and Hf (all FWHM-distinct significant found features
   outside every requested-source core, fitted jointly).  Linear, Quadratic, FlatStep, and LinearStep
   continua are re-fit for every hypothesis.  Conditional/quasi-AICc counts every continuum parameter,
   the one Hs scale, and every Hf amplitude.  Upstream found-peak means/widths are conditioned on rather
   than re-estimated, and weighted least-squares estimates are scored by measured Poisson deviance; this
   is deliberately not described as a fully free maximum-likelihood AICc.
2. Rejected provisional groups do not mint atoms and therefore cannot act as artificial bridges.  The
   rejection set is transactional: at least two supported groups must survive, otherwise all local
   rejections are rolled back.
3. For a remaining over-wide union, the whole component is scored at every atom-core-safe channel edge.
   The challenger is exactly two channel-disjoint peak-plus-continuum ROIs on the incumbent union's
   channels.  Its measured-data AICc includes the existing configured soft-width pressure; there is no
   hard splitter or new width constant.
4. The selected parent, two exact channel-aligned children, calibration, and atom ownership are carried
   as provenance.  They are captured before ordinary interference-edge shrinking, translated into the
   original calibration through channel coordinates, and applied one component at a time.  Each local
   solve must preserve recovered predicted anchors, significant clean-source anchors, and every
   incumbent observable requested-source peak.  A failed component rolls back without vetoing an
   independent accepted component or importing unrelated re-clustering geometry.
5. Diagnostics are finalized over the actual incumbent ROI range affected by the transaction.  Every
   overlapping H0/Hs/Hf, pair, and whole-component diagnostic is labeled as an accepted transaction or
   a rolled-back proposal; rollback reports `MergeInseparableWide` while retaining the proposed decision
   and reason as evidence.  Preliminary trace lines remain chronological proposals and are followed by
   the explicit final component disposition.

No blanket found-peak protection was retained.  The earlier fixed-Z external-attenuation candidate,
unconditionally broad geometry candidates, and broad anchor protection remain rejected.  Intermediate
geometry builds reached 61--78 full-corpus misses but introduced new U235/Ir192/Th232/Np and other
regressions; a partially guarded build still had 76 misses.  The root cause was applying whole-component
geometry before a collapse existed, allowing rejected local geometry to fall through to global
refinement, and coupling multiple component partitions into one all-or-nothing solve.  Source-collapse
gating plus exact per-component transactions removed those regressions.

### Am241_Sh mechanism

Final trace: `/private/tmp/shielded-source-recovery-roi-final13/traces/Am241_Sh.{out,err}`.

- Source-clean recovery retains 17 anchors and activity 108.2 MBq, within a factor of two of the
  126--131 MBq legacy reference rather than the original order-of-magnitude collapse.
- H0/Hs for 300.46 keV are 25.5/27.9 and for 312 keV are 40.2/42.6, so both unsupported groups are
  rejected.  The 322-, 332-, and 335-keV groups beat H0 and remain.
- The target whole component scores union AICc 726.2 versus partition AICc 723.8 with 9.6 width pressure.
  Its exact local transaction replaces the affected 293.5--348.1-keV incumbent with fit ROIs
  `[315.0,327.1]` and `[327.1,341.3]` keV.  The observable refit tightens these to approximately
  `[318.2,326.9]` and `[328.4,339.5]` keV.  There is no fitted 300/312-keV peak; 322.6 keV is isolated
  from the 332.3/335.3-keV cluster.
- The independent 93--133-keV proposal is rolled back because its local solve loses the real
  114.1-keV observable anchor.  It does not prevent the high-energy component from being accepted.
- Final output is 24 observable / 55 fitted peaks.  The eight new Am241_Sh misses from E006 are all
  recovered, including 322/376/383 keV through observable refit; the one remaining Am241_Sh definite
  miss predates this recovery.

### Final panel and corpus evidence

Stable 11-spectrum panel and gallery:
`/private/tmp/shielded-source-recovery-roi-final13/diagnostic/`.

- 11/11 Success, 3 definite misses, 0 spurious, 0 failure, scalar cost -150.155776563.
- Versus the retained root fix, the other ten spectra are byte-identical; Am241_Sh has the same
  miss/spurious result, four fewer private fitted peaks, one fewer observable peak, two more final ROIs,
  and a 14.64-keV / 7.22-FWHM reduction in the formerly dominant target union.  Scalar cost changes by
  +0.365 because the objective does not directly reward this visibly necessary geometry correction.
- Versus atom-safe E006: 11 -> 3 misses, 0 spurious unchanged.  Versus E000R: 4 -> 3 misses, 0 spurious
  unchanged, with scalar cost improving by 2.373.  Comparisons are `compare_e006_atom_safe/`,
  `compare_retained_root/`, and `compare_e000r/`.

Stable comparable 205-spectrum Detective-X/Baltimore/300-s corpus and complete gallery:
`/private/tmp/shielded-source-recovery-roi-final13/full_baltimore_300/`; comparison to the retained root
fix: `/private/tmp/shielded-source-recovery-roi-final13/full_compare/`.

- Definite misses 71 -> 70, spurious 7 -> 7, Success 200 -> 200, legitimate empties 5 -> 5, mechanical
  failures 0 -> 0, scalar cost -2878.16818763 -> -2881.92299705.
- No individual spectrum gains a definite miss or spurious peak.  Pu238_Sh loses one definite miss, so
  the retained geometry improves beyond Am241_Sh.  Named Am241_Unsh, Eu152_Sh/Unsh, Cs137_Unsh and the
  other ten panel sentinels retain their prior hard results.
- The exact-provenance/diagnostic corrections produce byte-identical panel and full result TSVs to the
  preceding accepted geometry build, while the full HTML gallery was regenerated from the corrected
  executable for user review.

### Isolation, review, and residual risk

The current R6/legacy Am241 observable N42 is byte-identical to the frozen baseline after normalizing
only `n42DocDateTime`.  The nondeterministically ordered private-fit N42 has identical structural
multisets of all 20 continua and all 305 peaks.  The four named non-collapse policy sentinels were also
byte-identical after timestamp normalization before the final provenance-only correction, and the final
panel/full TSVs remain byte-identical across that correction.  R2 admission and protected mixed/user
geometry continue to be exercised by the full source-fit suite.

The final adversarial re-review closed both medium findings: selected children are now captured before
unscored shrinking, and rollback/acceptance labels cover every related diagnostic over the actual
affected incumbent.  It found no remaining release-blocking code issue.  A low statistical caveat
remains: Hf is a conditional structural challenger because upstream peak locations/widths are reused;
the strong corpus gate, narrow activation condition, and final observable refit bound that risk.

Final verification on this exact build:

- focused `StatisticalDetailHelpers`: 22/22 passed;
- `PeakFitImproveReporterSelfTest`: 1/1 passed;
- full `test_fitPeaksForSources` with developer checks: 47/47 passed, `No errors detected`
  (`/private/tmp/shielded-source-recovery-roi-final13/full_tests.out` and `full_tests.err`);
- `git diff --check`: clean; nothing staged or committed.
