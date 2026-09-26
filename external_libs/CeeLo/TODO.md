# CeeLo — Open Items / TODO

Tracker for known shortcomings, deferred features, and fidelity/speed improvements.

**Process:** when an item is resolved, **delete it from this file in the same commit that fixes it**

Each item lists **Effect** (what's wrong / what it buys), **Implementation** (what it would take), and
**Shows up when** (the regime where it matters). `[severity]` is relative to the ~1%-vs-GEANT4 goal.
Items measured to be sub-statistical in the gated configs are kept for awareness, not because they
currently fail validation.

---

## Cascade / true-coincidence summing

### Rejected-branch pair links cannot all be honored `[high — decay data]`

**Effect.** Rejected level schemes fall back to flat member marginals plus directional pairwise
coincidence links. Those inputs are not jointly consistent: a whole-database census found 10,355
linked pairs whose implied joint exceeds `min(P(a),P(b))`, and 3,172 of 131,851 material pairs whose
two directional reconstructions disagree by more than 1%. FullRealization now projects those inputs
onto a deterministic one-parent Bayesian forest: every selected pair joint is Fréchet-bounded and
every member marginal is preserved, but incompatible links necessarily go unused. Conditional and
analytic estimators deliberately do not invent gamma-gamma sum-fed joints for rejected branches.
All three set per-peak `summing_model_complete=false` when an invalid primary or a materially possible
invalid-branch photon pair can affect that requested peak.

**Implementation.** Recover authoritative level topology or a consistent multivariate decay law for
the rejected records. The bounded forest is a safe realization, not a claim that discarded links are
false. Merely averaging directional conditionals or clamping every pair independently cannot close
this item because neither defines a coherent 3+-member realization. Analytic and Conditional should
enumerate the selected forest's pair joints only if the projected policy becomes part of their public
physics contract; until then the incomplete flag prevents silent use as an exact correction.

**Shows up.** Any peak marked `summing_model_complete=false`, especially rejected U-235, Pu-241,
Sb-125, and invalid branches whose two lower-energy photons can sum into a requested line.

### Infeasible ICC/intensity pairs can cost branches their level scheme `[medium — decay data]`

**Effect.** Some records pair a gamma intensity with an ICC for which
`I_gamma*(1+alpha) > 1`, i.e. the inferred transition would occur more than once per decay. The
adapter bounds the selected transition and rejects a graph whose reconciled feeding remains
infeasible. Counts in the original audit snapshot (1052 topology-bearing violations, 1848 expanded
occurrences, and 1897 distinct parent->child candidates -> 1542 valid / 355 rejected) are retained in
the study artifacts for provenance, but are **not current effective-tree counts**: they used the old
XML and an earlier acceptance policy. The final regenerated-XML test reports per-branch-instance,
disjoint audit categories: **4504 valid, 817 raw-only, 1782 partial-only, and 1328 both raw+partial**;
accepted graphs contain **323 E0 repairs and 801 intensity caps**. Those categories deliberately do
not collapse into the historical distinct-pair denominator; a final distinct parent->child recount
is still needed before quoting a replacement rejection percentage. Rejected branches fall back to
the pairwise-coincidence model described above.

Historically rejected examples with parent half-life > 1 y include **U-235 -> Th-231
(current raw_total_feed 2.37344), Th-229 -> Ra-225 (2.13),
Hf-172 -> Lu-172m (2.01), Pu-241 -> U-237 (1.96), Pa-231 -> Ac-227 (1.58), Sb-125 -> Te-125 (1.05)**,
Am-242m, Cf-249, Bk-247, Cm-243, Np-235, Ho-166m, Po-208, Sn-126 (two pairs), Es-252. Shorter-lived
but notable: Pa-234 -> U-234 (3.31, T½ 6.7 h) and Pa-234m -> U-234 (1.10, T½ 1.16 min). Counting
distinct parent->child pairs, a further 36 have T½ > 1 d and 41 T½ > 1 h. **Pu-241 matters for any
reactor-grade plutonium and
Sb-125 is a common fission product** — both silently switch estimator model. U-235's summing factors
change by 12–850x across the switch (measured in the dev-only ICC feasibility study; `u235_before.csv` vs
`u235_after.csv`), so no absolute magnitude derived from the level-path result for a rejected
nuclide is trustworthy. Nothing in the G4-anchored cascade gate
(`tests/data/geant4_reference/cascade_summing_multi.csv` — Am241/Ba133/Co57/Co60/Eu152/I125/Na22,
every one with entry_probability >= 0.9985 except Eu-152 -> Gd-152 at 0.99847) can detect a
regression here.

**Not a matching failure.** The G4 <-> SandiaDecay match is *correct* for the dominant offenders:
U-235's 19.59 keV gamma is correctly matched to Th-231 level 5 -> 4 (`z90.a231`, alpha = 114.7), and
its `cascade-match` audit remark reports no ambiguity. The original input evaluations disagreed —
the old SandiaDecay I_gamma = 0.61 with GEANT4's alpha = 114.7 implied 70.6 occurrences per decay.
The final regenerated XML carries the curated ENSDF override instead.

**What is already done.** Emission is bounded rather than guessed: a member's gamma intensity is
capped at `1/(1+alpha)`, the hard ceiling that holds whichever tabulated number is wrong, so no
branch — accepted or rejected — can emit a line above what its transition supports. That degrades
gracefully where a semantic reinterpretation does not: I-125's saturating Te-125 transition
(I = 0.0668, alpha = 14.08 -> 1.0073) loses 0.7% while U-235's 19.59 keV loses a factor of 70,
removing a 21x line excess and the false sum peaks it fed onto 205.31 and 221.40 keV.

**What remains.** The final U-235 override puts the 19.59-keV photon at the ENSDF central value,
**0.00583 per decay**; the earlier 8.64e-3 / 1.5x statement no longer applies. For records without a
curated override, the generic physical bound is still a ceiling rather than an estimate. Reading a
tabulated intensity as the TRANSITION rate is right for some records (Ag-104m's 6.9 keV isomeric
transition lists intensity exactly 1.0 against alpha = 1e8) and catastrophic for others. The
historical snapshot had 435 of 1052 violations within 2% of the bound, where the two readings differ
by 15x. Nothing but violation size separates them, so improving the uncurated cases needs a
per-transition ENSDF cross-check of which quantity the intensity is, not another heuristic. Census tool:
the dev-only `icc_feasibility_sweep` probe.

**Shows up.** U/Pu chain nuclides with strongly-converted low-energy transitions; harmless for the
calibration set (Co-60, Ba-133, Am-241, Eu-152, Co-57, Na-22 are all unaffected).

### Internal-pair formation for verified high-energy E0 transitions `[low — physics]`

**Effect.** Authoritative primary-evaluation provenance now marks verified pure-E0 transitions
independently of GEANT4's advisory multipole code. Every verified E0 gets `p_gamma = 0`; below the
exact 2m_e threshold its unpersisted shell remainder is approximated by one shell-unresolved
conversion electron. Above threshold the remainder may be INTERNAL PAIR FORMATION. It is retained as
an explicit `p_unmodeled` outcome so the level path advances without inventing a photon, conversion
electron, vacancy, or deposit. N-16's verified 6048.2 keV 0+ -> 0+ transition is therefore no longer
emitted as an impossible photon, but its real pair deposit is still absent.

**Implementation.** Emit a back-to-back 511 pair plus the (E - 1022) keV kinetic energy for the
internal-pair branch, with the pair/conversion split from data where available. Do not promote
GEANT4 code 1 alone to an E0 assertion: Ac-232's observed 373.3 keV and Rb-78's observed 189.8 keV
photons have blank/ambiguous ENSDF M/CC assignments despite code 1 and must remain photons.

**Shows up.** Verified E0 transitions above 1022 keV.

### Na-22 511 keV drifted from GEANT4 after the beta+ model change `[medium — accuracy]`

**Effect.** The cascade-summing snapshot `our_cascade_multi.csv` (now regenerated, not
committed) was last written at `8ca66d8` (2026-07-07). Re-running its producer (`examples/cascade_observables`) at the same 8M/4M statistics on
the current main tree gives **Na-22 P511 +3.24% (z = +23.8)** and sum1785 +3.81% (z = +6.0);
everything else is within |z| <= 3.2. Engine/G4 agreement for Na-22 P511 moved 1.0012 -> 1.0340.
Bisects to **d3e0ee3** (2026-07-28), which replaced `emit_annihilation` with a full beta+ spectrum +
positron-range + in-flight-annihilation model. The G4 reference is unchanged over that span, so this
is entirely engine-side. Still inside the gate's band, so nothing failed.

**Why it went unnoticed.** `profiling/compare_validation.py` is the only consumer of that CSV and is
**not wired into ctest** — the cascade dashboard is manual, so a 24-sigma drift sat for three weeks.
Worth adding the cascade block to ctest.

**Implementation.** Triage the positron model against G4 for Na-22, then regenerate the reference
(`cascade_observables`, ~85 s) and record the producing commit in its header. Do not regenerate first
— that would bake the drift in as the new expectation.

### Sum-peak-fed / x-ray-fed summing-IN in the Conditional estimator `[low]`
- **Effect:** `CascadeMethod::Conditional` now captures gamma-gamma sum-peak-fed summing-IN (pairs a+b
  in the window where neither is the peak gamma, e.g. Co-57 122.06+14.41→136.47) and coincident vacancy
  K/L x-ray summing-OUT (via `cascade_sum_pair_channels` / `cascade_level_vacancies`, in
  `cascade_peak_thread` — the default method is no longer γ-only for x-rays). Residual vs FullRealization
  at contact: **≈1–2% on Co-57 136 (1.18 vs 1.20) / Ba-133 356 (0.517 vs 0.494)** from (a) x-ray-FED
  feeding not enumerated — γ + coincident K x-rays adding into the window (Co-57 122 + two Fe Kα, where
  14.41 converts), and (b) the partner-independence approximation for 3+-gamma cascades. ε²-small at far
  geometry (Co-57 136 @10 cm is 1.020 vs 1.022, within stats).
- **Implementation:** enumerate γ + x-ray-pair (and triple-γ) sum-fed channels in
  `cascade_sum_pair_channels`; and/or sample the peak's partners jointly (level-path realization) rather
  than independently. Or route contact-geometry x-ray-heavy summers (Ba-133, Am-241) to FullRealization.
- **Shows up when:** high-Z low-energy X-ray summers or sum-peak-fed peaks at contact with the default method.

### Analytic cascade-summing (`compute_cascade_analytic`) remaining approximations `[low]`
- **Effect:** the fully-analytic path (`src/cascade/AnalyticCascade.{h,cpp,_imp.hpp}`, the InterSpec DRF entry point)
  matches FullRealization to ≲1% at all geometries incl. contact for the mainstream γ-γ / EC-x-ray / β⁺
  cases (analytic-vs-FR MC @2 cm: Co-57 122/136 z≤0.4, Ba-133 356/302 z≤0.97, Na-22 1274 z=0.68,
  Am-241 59.5 z=0.43, Eu-152 344/1408 z≤0.74, Co-60 z<0.5; far z<0.15). It uses the exact level-scheme
  survival DP (Eq. 10′) for summing-OUT and a unified γ/x-ray/511 occurrence enumerator (pairs + triples,
  no MC triple-subtraction, same-vacancy exclusivity) for summing-IN — each sum-fed channel additionally
  carries a **coincident-survival factor** (the summing-OUT of the sum-fed pair by everything else
  co-emitted). Back-to-back 511 pairs are summed-out with **2·ε_tot** (not the independent 2ε−ε²): for a
  single-sided detector at most one of the pair can deposit (Na-22 1274 z 1.7→0.68). **Overlapped windows**
  (multiple emitted lines in one peak window, e.g. Ba-133 79.6+81) are handled with a fitted-peak-area SF
  (emission-weighted mean over the window's lines). Known residuals:
  - **Low-energy partial-deposit / escape-peak recovery (the ε-only summing-out limit).** `Π(1−p·ε_tot)`
    treats the peak as a delta at the FEP and the coincident photon as binary remove/keep; it misses events
    where the primary's *escape peak* + the coincident x-ray sum back INTO the window. Worst case measured:
    **I-125 35.49 @2 cm, −4% (analytic over-removes)** — 35.49 sits just above the NaI iodine K-edge (33.17),
    so it K-escapes to ~6.9 keV and the coincident Te Kα (27.5 ≈ the 28.6 keV iodine escape energy) sums it
    back in. Also **Ba-133 81 @2 cm (−6%)** — low-energy, fed by the whole cascade above it + EC Cs K x-rays;
    milder on Eu-152 121.78 (−1.5%). This is the accepted limitation of the standard total-efficiency summing
    method (GESPECOR et al.); it is **NOT fixable from ε_FEP/ε_tot alone** — the provider would have to expose
    the deposit spectrum (escape-peak energy + fraction). ~0 by 20 cm. I-125 35.49 and Ba-133 81 @2 cm are
    report-only in the gate for this reason.
  - **Kα→L secondary vacancy** not modeled (a K vacancy relaxing via Kα leaves an L vacancy that also
    radiates; FR emits it). Deliberately dropped: for summing-OUT it is idempotent (the K x-ray already
    removed the event) and for summing-IN ε²-small — it affects only the sum *continuum*, not peak SF
    (Ba-133 Cs Kα→L agrees <1% without it).
  - **W(θ) angular correction** uses the collinear limit **g = W(0) = 1 + a2 + a4** for correlated
    summing-IN pairs (default ON, geometry-free): coincident FEP detection requires both photons to point
    into the detector, so their mutual angle is small and W(θ_ab)≈W(0) (the FEP acceptance is on-axis-peaked;
    a uniform-cap average under-weights the small angles — an earlier disk-subtense model measurably *hurt*
    at contact). Recovers Co-57 122+14→136 to ~0.2% (vs ~1.9% at g=1). Summing-OUT partners get g=1 (they
    need only ANY deposit; effect negligible — Co-60 <0.5%). Residual: W(0) is the collinear upper bound;
    the true FEP-weighted average is slightly below it, so strongly-correlated sum-in at very close geometry
    (large cone) may be marginally over-corrected — bounded by the a2/a4 magnitude (<W(0)−1, i.e. ≲ a few %).
  - **Triple-fed summing-IN** is enumerated only for **gamma lines** on distinct transitions (order 3),
    sorted + energy-pruned; x-ray-line and EC/annihilation triples and 4+-member sums are omitted. This is
    ε³-and-branch-small AND a perf necessity: enumerating over the ~25 expanded vacancy x-ray lines per
    converting transition is O(n³) and was ~27 s on Eu-152 / minutes on Am-241 before the restriction (now
    ~0.1 s / ~1 s, dominated by the SandiaDecay `build_cascades` call, not the analytic).
    Gamma products retained as topology-free categorical residuals participate in analytic and Conditional
    **pairs**, but are not yet admitted to this optional analytic triple enumerator; a peak fed only by a
    residual-containing three-gamma sum therefore misses that epsilon-cubed term.
  - **IC-electron deposition** is inherited from whatever the ε provider models (the analytic path only
    combines ε_FEP/ε_tot); it does not add the IC-electron summing that `compute_cascade` now models via
    `enable_ic_electrons` (see DESIGN.md → Purpose, cascade summing corrections).
  - **Non-gamma primary peaks.** A 511 primary uses a marginal-product C_out fallback (partner 511 +
    whole-scheme survival) — less rigorous than the gamma survival DP, and FullRealization's locator can't
    score a 511 peak, so it is validated only indirectly (via the coincident-511 effect on Na-22 1274). A
    peak requested at a pure **vacancy x-ray** energy (e.g. an Am-241 Np L line) returns found=false: the
    x-rays are generated from vacancies, not carried as `members`, so there is no line to locate — x-ray-line
    peaks as a *primary* are unsupported (they still contribute to other peaks' summing as coincident lines).
- **Implementation:** the escape-recovery term needs a deposit-spectrum-aware provider (escape energy +
  fraction); a real angular joint-acceptance from the DRF; extend the enumerator to EC/annih triples if a
  case demands it.
- **Shows up when:** very-low-E lines (≲60 keV) just above the detector K-edge with a coincident x-ray near
  the escape energy, at contact; or sub-0.5% agreement for high-Z x-ray-rich EC nuclides at contact.

### Coincident L x-ray sum peaks slightly low — upper-level feeding `[low]`
- **Effect:** with the per-subshell L + Coster-Kronig vacancy model, L x-ray *singles* match GEANT4 to
  ~0.3% (Am-241/Np, Pb-203/Tl, source-region-matched geometry), but the *coincident* L sum peaks run a
  little low: Am-241 59.5+L ≈ 0.80, L+L ≈ 0.96 vs G4. Since the singles (dominated by the anti-coincident
  59.5-IC vacancy) match, the deficit points to the upper-level cascade feeding *into* the 59.5 level
  being slightly under-weighted → too few 59.5-coincident upper-cascade L x-rays.
- **Implementation:** check the flow-conservation feeding in the level-path builder
  (`SandiaDecayCascade.cpp`, `feeding = out_flow − in_flow`) against the known Am-241 α-branching to the
  Np 102.96 / 158 keV levels; the discrepancy is likely in how direct (α/β/EC) level feeding is inferred
  vs gamma in-flow when intensities are imperfect.
- **Shows up when:** sum-*peak magnitudes* for high-Z L-x-ray emitters at close geometry. Negligible for
  peak-efficiency summing corrections (59.5+L ≈ 1.2% of the 59.5 peak → <0.25% on its efficiency).

### IC-electron summing — residuals & auto-enable follow-ups `[low]`
- **Effect:** `compute_cascade`'s `enable_ic_electrons` (default OFF; K/L conversion + K-Auger electron
  deposition, air/distance-gated, GEANT4-validated to <1% at 2 cm — see DESIGN.md → Purpose, cascade summing
  corrections, and the dev-only IC-electron GEANT4 cross-check). Two known gaps:
  - **(a) M+/outer-shell and sub-15 keV conversion electrons are not deposited.** These carry ~the full
    transition energy but are dropped, slightly under-counting summing-OUT on the IC-heavy low-E peaks at
    **bare-source contact / 4π** (Ba-133 P356 ≈ +14% vs G4 there); ≈0 by 2 cm as air / real source
    material absorbs them. A crude full-energy M+ deposit was **tried and reverted** — it mis-shifts the
    converting transition's OWN low-E peak and over-corrects at 2 cm (1.00→0.96). Needs proper M-shell
    binding/energy + partial-deposit (backscatter/escape) modeling, not a flat full-energy dump.
  - **(b) Not auto-enabled by proximity.** The flag is default-OFF and costs ~5–8% when ON regardless of
    distance (the per-conversion isotropic-sample + detector ray-trace + air-range CSDA lookup runs even
    when the electron provably can't reach). An `Auto` mode would give correct-by-default close-source
    physics and skip the cost far away.
- **Implementation:** (a) add M-shell binding data + a partial-deposit factor (or a tabulated conversion-
  electron energy/deposit model). (b) replace the `bool enable_ic_electrons` with a tri-state
  `Off/On/Auto` (default `Auto`), enable within ~a few cm (where ΔSF > ~1%) — but re-tune the cascade
  test bands, since turning IC on at 2 cm shifts the (currently OFF-tuned) gated peak areas toward G4.
- **Shows up when:** (a) bare / thin-source contact / near-4π low-E x-ray-rich summers; (b) close-source
  callers wanting IC summing without setting the flag.

### L1/L2 conversion electrons use the L3 binding energy `[low]`
- **Effect:** both `build_radioactive_emissions()` and `compute_cascade`'s IC-electron deposition take
  the L-shell binding from `LFluorescenceData::l3_edge_keV` for L1, L2 and L3 conversion alike, so L1/L2
  conversion electrons come out too energetic by B(L1)−B(L3) or B(L2)−B(L3): up to ~4.8 keV for Np
  (L1 22.4 vs L3 17.6 keV), ~2.8 keV for U. Found Sep 2026 while fixing the daughter-Z gate.
- **Implementation:** carry per-subshell L binding energies in `LFluorescenceData` (the EADL generator
  already reads them) and use the matching one per `ConversionShell`.
- **Shows up when:** IC-electron summing or source-escape studies of strongly L-converted heavy-element
  transitions (Am-241, U/Pu chains); the energy shift is a few keV on 10–60 keV electrons.

### β⁺ 511 annihilation modeling — residuals `[low]`
- **Status:** the cascade β⁺ path now RANGES the positron (commit 81121e3): KE sampled from the allowed
  β⁺ spectrum (endpoint carried on the Annih511 member), (1) in-flight annihilation → two non-511 γ,
  (2) escape from the source → no clean source 511 pair; a contained positron annihilates at rest near
  the vertex. This took Na-22 P511 from +12% to +7% vs GEANT4 (matched geometry).
- **Remaining residuals:**
  - **Dominant (≈+7% on Na-22 P511):** the source walk UNDER-predicts positron escape, so too many clean
    source 511 pairs survive. This is the Molière-walk accuracy problem (see **B2/B3**) — NOT the
    annihilation model — compounded by the Fe-calibrated albedo gate possibly mis-applied to low-Z Al.
  - **(3) positron-specific stopping** not applied — the walk uses the *electron* term (Møller F⁻); the
    positron term (Bhabha F⁺) is ~80 lines, default-off (`bool is_positron` threaded through
    `walk_in_source_geometry`), **<1 %** effect, no perf cost (the walk already runs). Low leverage.
  - **No ortho-positronium 3γ** (quenched to 2γ in dense media; minor); **Conditional estimator** still
    point-annihilates (the ranging is in FullRealization only).
- **Shows up when:** β⁺ emitters (Na-22, F-18, …) in close/contact geometry — the 511 / 1022 / γ+511 peaks.

### Fully-independent Conditional-summing-factor-vs-GEANT4 gate `[low — validation]`
- **Effect:** the cascade-summing gate is now surfaced through `profiling/compare_validation.py` and the
  C++ ctest, both driven from one reference CSV (`tests/data/geant4_reference/cascade_summing_multi.csv`),
  and includes `estimator=cond` rows. But the Conditional rows are gated as a per-decay area
  `A_FR·k_cond/k_full` — i.e. anchored to the FullRealization per-decay normalization, with only the
  Conditional/FR summing-factor RATIO independent. A truly independent Conditional-summing-*factor*-vs-G4
  check would need a G4 "no-summing" baseline (the peak's per-decay area with coincidence summing removed)
  to divide the correlated area by; the harness has no such mode today.
- **Implementation:** a "no-summing" G4 baseline via **time-staggering** — one correlated run yields BOTH
  `A_corr` (with-summing) and `A_nosum` (no-summing), so `k_G4 = A_corr / A_nosum` and the engine
  Conditional `summing_factor` compares directly to `k_G4` (estimator-independent; the G4-vs-SandiaDecay
  emission-intensity difference cancels in the ratio).
  - **Idea:** the global time is used as a per-decay-gamma *lineage tag*. Offset each true decay gamma to
    its own ms-wide time window; its EM shower (recoil e⁻, fluorescence/brems x-rays) inherits the offset
    automatically (children take the parent's global time at creation), so each gamma's FULL deposit —
    hence its full-energy peak — stays together while decay-*sibling* gammas are pushed a full ms apart.
    Shifting the start time changes no physics (no time-dependent cross sections); it only relabels.
  - **Harness (`tools/geant4_validation/src/`):** in a **`G4UserTrackingAction::PreUserTrackingAction`**
    (cleaner than a stacking-action `const_cast`), stagger ONLY the genuine decay gammas —
    gate on the creator process, NOT "ParentID>0 && isGamma":
    ```cpp
    auto* cp = track->GetCreatorProcess();
    if (track->GetDefinition() == G4Gamma::Definition() &&
        cp && cp->GetProcessName() == "RadioactiveDecay")
        const_cast<G4Track*>(track)->SetGlobalTime(
            track->GetGlobalTime() + track->GetTrackID() * 1.0 * CLHEP::ms);
    ```
    Gating on `RadioactiveDecay` is CRITICAL: offsetting *every* secondary gamma would fling a primary's
    own fluorescence/brems photons into a different window and corrupt that primary's FEP.
  - **Scoring:** the time offset does nothing alone — change `EventAction`/`SteppingAction` to bin
    `(t_prestep, edep)` into ms-wide windows and emit ONE histogram entry per occupied window (= the
    no-summing spectrum). Summing all windows in the event reproduces the with-summing spectrum, so both
    fall out of the single run. Gate behind a `--no-summing` flag; enable per-config in the CSV
    (`estimator=cond` rows carry a directly-comparable `k_G4`).
  - **Caveats:** β⁺ / annihilation pairs (Na-22 511+511 share a positron ancestor, not a decay-gamma
    creator process → still co-window) and coincident x-ray lines need an extra lineage rule; the pure
    γ-γ / γ-xray summers (Co-57/Co-60/Ba-133) are clean. Confirm no global-time tracking cutoff is active.
    Nuclides beyond the current set (Cs-137 662-conv, Tc-99m) bake new alcyl rows into the CSV the same way.
- **Shows up when:** wanting an estimator-independent Conditional check rather than the FR-anchored ratio.

---

## Electron / CSDA physics

### B1 — in-crystal positron annihilation — residual approximations `[low]`
- **Done:** the PP positron is now CSDA-tracked and annihilates at its endpoint, not the vertex
  (`deposited_in_scoring` returns `stop_position`); it can annihilate **in flight** (per-substep Heitler
  rate → a photon pair of 2·mₑc² + residual KE, `annihilated_in_flight`); it walks with the
  positron-specific **Berger-Seltzer F⁺** stopping (`range_table_pos_`, `is_positron` flag); and the
  e⁺/e⁻ get the **Koch-Motz / G4ModifiedTsai** pair opening angle. Measured (2614 keV NaI, 8M A/B):
  in-flight annihilation drops the single/double-escape peaks ~3% — taking them from +2–3% (z≈+7) ABOVE
  GEANT4 onto it (double-escape −0.5%/z=−1.2, single-escape −1.3%/z=−4.5); the endpoint/stopping/angle
  pieces are sub-mm and statistically invisible in FEP/total (|z| ≤ 1.3, incl. thin CZT). FEP/total
  unchanged throughout.
- **Remaining (minor):** the in-flight pair uses the isotropic energy-split approximation (same as the
  cascade β⁺ path), not the exact 2-body boost kinematics (correlated angle/energy) — a candidate for the
  small single-escape overshoot; no ortho-positronium 3γ (quenched to 2γ in dense media anyway);
  `positron_annih_xsec_per_electron` is duplicated in `ElectronCsda.cpp` and `EfficiencyCalculator.cpp`
  and should be unified into one physics home.
- **Shows up when:** escape-peak / annihilation-continuum shape studies at E > 1.5 MeV, thin/edge crystals.

### B2/B3 — retire the empirical skin-escape gate via a Mott/G-S MSC `[low — robustness; research]`
- **Status:** the gate is **regime-aware** (`source_escape_survival_exit`, exit-state + light-Z floor +
  exit-energy window) — replaced the old birth-energy gate + positron-only flag + Fe-only `kAlbScale`,
  validated Z = 5→82 to ~1–2% vs GEANT4. Three escape models now live behind a **compile-time switch**
  (`-DCEELO_SOURCE_ESCAPE_MODEL=gate|tails|gs`, default `gate`; env A/B knobs
  `MCDET_NO_{STRAGGLE,B2,TAILS}`). Full per-config comparison: DESIGN.md → Known Limitations →
  "Source-electron skin escape".
  **Verdict (mean |err| @662): gate 0.56% < tails 0.62% < gs 0.92%** — none of the first-principles
  variants beats the gate; `gate` stays the default.
- **Why tails/gs fall short — BOX-specific:** the per-step Gaussian-core + screened-Rutherford tail
  **over-diffuses box geometries** (a large-angle scatter is as likely to redirect toward a wall as away).
  This is geometry-specific: on sphere/shell geometries (G-B/G-C, cfg 11) all three models agree to ~±1%
  — the over-escape only appears for multi-wall boxes. A double-counting
  variance subtraction halved the `tails` excess (cfg 12 @662 +2.10→+1.07% vs gate +0.46%). `gs` sets the
  soft core from the screened-Rutherford soft transport moment (instead of Highland) but **under-scatters
  ~30%** (the screened cross-section lacks the **Mott** spin + nuclear-form-factor corrections), making
  the mid-Z boxes worse (+1.65%). A uniform G₁ scale-up just reproduces `tails`.
- **Path to actually retire the gate:** replace the screened-Rutherford cross-section with **Mott** (or a
  tabulated **Goudsmit–Saunderson** angular distribution sampled by inversion — no Gaussian+tail
  double-count). Runtime is ~free (the walk is <10% of total; a table lookup ≈ current cost); the cost is
  developer time (hundreds of lines + Mott data) + re-validating ~1% across Z without the gate.
- **Shows up when:** shielded/extended box sources (cfg 8/11/12/20-24); the gate remains the default.

### W/Pb high-Z box common offset `[low — investigate]`
- **Effect:** the new W (cfg 21) and Pb (cfg 20) wall-box configs show a ~+2–3.7% total over GEANT4 at
  2614 keV in **both** the gate AND principled modes (a *common* offset, so NOT the electron-escape method).
  Suspected: source-shield secondary channels (PP 511s / shield brems) or photon transport in very high-Z
  walls, exposed only now that cfg 20/21 exist. These are out-of-sample configs (no historical reference).
- **Implementation:** decompose with `--no-source-brems` / `--no-source-electrons` and `--entry-diag`;
  generate higher-stats G4; check the W/Pb fluorescence + PP-annihilation channels.
- **Shows up when:** Z ≳ 74 wall boxes at high energy (cfg 20/21).

### Electron tables for Np–Cf are uranium's `[low — coverage]`
- **Effect:** ESTAR stopping, ICRU-49 I and the NIST EPQ Seltzer–Berger tables stop at Z = 92, so Z 93–98
  reuse uranium's (`electron_table_z()`). Measured against ESTAR 93–98 (DESIGN.md → Known Limitations):
  collision stopping within −2.1 .. +2.3%, radiative −1.9% (Np) .. −8.9% (Cf), Pu −2.6 .. −4.5%;
  CSDA range −1.7 .. +2.9% above 50 keV. The effect on efficiencies has not been measured.
- **Implementation:** add ESTAR 93–98 as a production source (the NIST text interface serves them; they
  were fetched only to measure this approximation) and a bremsstrahlung source above Z = 92 (EPQ ends at
  pdebr92; EEDL2023 MF=26 MT=527 covers all Z but would be a new lineage). A cheaper step for the
  radiative part alone: rescale uranium's by Z(Z+1)/A, which brings the measured radiative error to
  −2.5 .. +3.9% (Pu −0.8 .. +1.3%).
- **Shows up when:** thick Pu/Am/Cf sources or shields at MeV energies, where bremsstrahlung and
  electron escape from the actinide contribute to the total efficiency.

### CZT electron escape `[resolved Aug 2026 — small residual]`
- **Effect:** was FEP −3 to −7% for thin crystals >1 MeV. The Aug 2026 crystal-walk fixes
  (path-consistent Highland + per-step Bohr straggling + step-budget guard) brought cfg 5 to
  +0.66% ± 1.11 @800, −0.31% ± 1.31 @1000, −1.31% ± 1.81 @1500 (all |z| < 0.7) — straggling was
  the missing physics for thin-crystal escape, not a marginal help.
  See `studies/high_e_fep/FINDINGS.md`.
- **Shows up when:** any residual is now below the ~1.3-1.8% MC+ref precision of the gate rows.

---

## Photon transport physics

### 10–20 MeV transport is not validated `[low — coverage]`
- **Effect:** the photon tables reach 20 MeV (Sep 2026) so that attenuation is available for high-energy
  reaction lines, but no GEANT4 comparison covers transport above ~3 MeV. Pair production splits the
  kinetic energy evenly between the leptons (spectral shape, not attenuation), and photonuclear
  absorption is absent — as in EPDL and XCOM — although at the giant resonance (12–15 MeV for heavy
  nuclei) it adds a few percent to the attenuation of Pb/U (≈0.64 b for Pb-208 at 13.4 MeV, against
  18.8 b atomic).
- **Implementation:** a GEANT4 gate row at 10–15 MeV for a bare and a Pb-shielded detector; sampled
  pair-energy sharing; per-element photonuclear cross sections if they ever matter (parked).
- **Shows up when:** reaction gammas above ~10 MeV, e.g. N-14(n,γ) 10.8 MeV.

### cfg 8 Marinelli FEP −2.1% at 59 keV (near-peak scatter) `[low — residual]`
- **Effect:** the *only* genuine cfg-8 residual after the June 26 2026 reference regeneration. The old
  "low-E Marinelli deficit" (FEP −4 to −7%) and the "−2.5% albedo-recovery deficit @662" were **stale-
  reference artifacts** (the March reference predated the June toolchain; resolved — reference regenerated
  at 32M on the current geometry, all six cfg-8 energies now gate with no SKIP). Against the fresh
  reference, MC agrees to ≤0.7% total / ≤1.1% FEP for 100–2614 keV; the residual is **FEP −2.1% at 59 keV**
  (z≈5). The G4-harness u/s decomposition localizes it: unscattered FEP `eps_u` matches G4 to 0.01–0.3% at
  every energy (geometry/self-attenuation correct), but the small *scattered* FEP stream `eps_s` is ~15%
  low — MC under-places forward Rayleigh / small-angle-Compton photons in the ±1.5 keV window. Shrinks to
  ≈0 by 200 keV.
- **Implementation:** improve near-peak forward-scatter recovery (Rayleigh/small-angle-Compton landing in
  the FEP window) at low E; possibly the same crystal-escape albedo-recovery model that would also catch
  the small residual real albedo effect (G4 recovers crystal-escaped e⁻/γ backscattering off the
  surrounding water). Low leverage — the effect is <2.1% and only at the lowest energy.
- **Cross-check (2026-09-04, bench cfg 27/28: GEM35-70 + point source + 0.5 cm Fe shell at 2/10 cm,
  60/88/122 keV, 0.75 keV window on both sides):** behind DEEP iron the absolute CeeLo-vs-G4 FEP gap is
  dominated by the documented EPICS2023-vs-EPICS2014 iron photoelectric difference (~2% in mu, 8% through
  4.6 mfp; CeeLo = NIST), so compare the u/s FEP SPLIT there, not the absolute: the scattered share of
  the 60 keV peak is 30.0% (CeeLo) vs 30.6% (G4), i.e. the forward-scatter placement agrees to ~2% of
  that stream in a Rayleigh-dominated case.  The analytic side of the same physics is
  `rayleigh_deflection_loss_fraction` (see DESIGN.md, FEP window section).
- **Shows up when:** Marinelli / wrap-around extended sources at ~59 keV (FEP only; totals agree to ≤0.7%).

### LaBr3 high-energy FEP residual `[low — residual]`
- **Effect:** after the Aug 2026 crystal-walk fixes (which closed the old "high-energy FEP family";
  see `studies/high_e_fep/FINDINGS.md`), cfg 3 (LaBr3 2"x2") FEP sits +0.63% (z +4.0) @2000 and
  +0.49% (z +2.6) @3000 at 0.16% MC precision — the only config-specific structure left ≥2 MeV
  (all NaI configs |z| ≤ 1.7; pooled ≥2 MeV +0.24%).
- **Implementation:** candidate causes: Bohr straggling variance slightly large for LaBr3
  (no shell/density corrections), or a small cfg-3 reference offset. Well inside the 1.5% tolerance.
- **Shows up when:** LaBr3 ≥1 MeV photopeaks, at the +0.5% level.

### Pb-attenuator total-efficiency excess `[low — residual]`
- **Effect:** cfg 7 (2 mm Pb attenuator) TOTAL efficiency +0.25 to +0.27% (z ≈ 3) at 2000/3000 keV
  at 0.08% MC precision (pre-fix +0.41/+0.63%, z up to 7.6); FEP agrees. Pb-specific, grows with E.
- **Implementation:** untriaged; plausibly attenuator pair-production/secondary treatment.
- **Shows up when:** high-Z attenuators ≥2 MeV, totals only.

### L/M-shell fluorescence X-ray escape not modeled `[low]`
- **Effect:** only K-shell fluorescence is transported; L X-rays (e.g. Pb Lα ≈10.5 keV) are deposited
  locally → over-predicts local deposit and misses L-escape spectral features near surfaces.
- **Implementation:** add L-line emission using the direct EPICS2023 EADL
  fixtures in `src/cross_sections/relaxation_epics_data.cpp`; the remaining gap
  is transport/use of those lines, not source data.
- **Shows up when:** thin high-Z detectors (CZT) or PE in a high-Z attenuator skin; low-energy spectra.

### Source-material fluorescence not emitted (high-Z self-attenuating sources) `[medium — coverage]`
- **Effect:** photoelectric absorption **in the source material / source shields** emits the
  photoelectron (for brems) but **no characteristic K X-ray** (`SourceGeometry.cpp`,
  `transport_source_photon_impl`, `InteractionType::Photoelectric` returns `survived=false`). So K
  fluorescence that escapes the source and reaches the detector is missing from **total** efficiency.
  Quantified June 2026 against GEANT4 for a self-attenuating **Thorium** sphere (3"×3" NaI): MC total is
  **−22% (238.6 keV), −4.1% (583 keV), −2.6% (911 keV)** below G4 — the deficit tracks the PE fraction
  (the G4 histogram shows the missing counts as a Th Kα line at 90/93 keV). **FEP is unaffected** (the
  escaped X-ray cannot land in the primary photopeak) — FEP agrees to |z|≤1 across the same runs — and
  **low-Z trace sources are unaffected** (soil sphere G-B/G-C totals agree to ≲1.5%). Orthogonal to
  geometry: it applies to any high-Z source shape, not just the sphere.
- **Implementation:** mirror the crystal-side K-fluorescence path — after PE in source material/shield,
  sample a K-shell vacancy (energy-dependent K fraction) and emit one fluorescence line as a
  `SourceSecondaryPhoton`, transported through the remaining source geometry like the brems secondaries
  (the K-line data is already in `relaxation_epics_data.cpp`). Gate it so low-Z (<10 keV) lines, which deposit
  locally anyway, add no cost; verify cfg 8/11/12 totals are unchanged.
- **Shows up when:** self-attenuating or thickly-shielded **high-Z** sources (U/Th/Pu metal, Pb-encased)
  at sub-MeV energies where PE dominates.

### cfg 7 100 keV Pb K-edge `[G4-side, excluded from gate]`
- **Effect:** ~−25–31% total at 100 keV — localized to G4's in-tracking μ_total(Pb) being smeared across
  the K-edge (G4-side, not an MC bug). Excluded from `compare_validation.py`.
- **Implementation:** none on our side; revisit if G4 fixes the edge interpolation. Separately verify our
  Pb K-line branching ratios (MC Kα/Kβ 0.67 vs G4 0.82) against EADL.
- **Shows up when:** just above high-Z K-edges behind thick high-Z attenuators.

---

## Cross-section data / fidelity

### Actinide atomic weights are one conventional mass per element `[low — coverage]`
- **Effect:** above uranium there is no standard atomic weight, so CeeLo uses xraylib's per-element
  masses: Np 237, Pu 239.1, Am 243, Cm 247, Bk 249, Cf 251. Pu and Bk are the isotope actually in
  hand (Pu-239, Bk-249; confirmed Sep 2026); Am, Cm and Cf are the longest-lived isotope. A `Material`
  turns element mass fractions into atom densities with these masses, so a material that is really
  Am-241 gets 0.8% too few atoms per gram and a μ 0.8% low (Cm-244: 1.2%; Cf-252: 0.4%).
- **Implementation:** let `Material` take a per-component atomic mass (or a nuclide composition) so a
  caller that knows the isotopic mix — InterSpec does — sets it directly; keep xraylib's values as the
  default.
- **Shows up when:** per-gram quantities for Am/Cm/Cf sources or shields, e.g. the self-attenuation of
  an AmO₂ source or activity from mass.

### Atomic weights through uranium are xraylib's two-decimal table `[low — fidelity]`
- **Effect:** `g_atomic_weights` holds xraylib 4.2.1's values, rounded to two decimals and partly out
  of date: H 1.01 (IUPAC 1.008, +0.20%), Au 197.2 (+0.12%), Ti 47.9 (+0.07%), Mg 24.32 (+0.06%),
  Ge 72.59 (−0.06%), Na 23 (+0.04%), Al 26.97 (−0.04%), and Tm 167.27, which is erbium's value
  (Tm 168.93, −0.98%). `Material` turns mass fractions into atom densities with them, so each
  element's atom density is off by the same fraction in the other direction: hydrogen in water 0.2%
  low (electron density −0.04%), germanium in HPGe 0.06% high. The GDML export has its own IUPAC
  masses, so CeeLo and GEANT4 already differ by these amounts — far below the gates' precision,
  except for thulium.
- **Implementation:** take the IUPAC 2021 standard (abridged) atomic weights for Z ≤ 92 as a locked
  source in `generate_element_support.py`, keeping the actinide conventions above; then re-run the
  migration validators and spot-check the GEANT4 gates. InterSpec's formula-to-mass-fraction
  conversion (`CeeLoUtils`) uses SandiaDecay's natural masses through uranium for this reason.
- **Shows up when:** hydrogenous and HPGe materials at the 0.05–0.2% level; anything containing
  thulium at 1%.

### Rayleigh anomalous dispersion (f′, f″) not included `[low]`
- **Effect:** the coherent cross-section uses non-dispersive F(x,Z) → overestimates by several-% to
  tens-% in narrow windows just above high-Z absorption edges.
- **Implementation:** add f′(E)+if″(E) (xraylib Fi/Fii or RTAB tables) to the Rayleigh XS / angular model.
  Check first (Sep 2026): against XCOM, which omits anomalous scattering, EPICS2023's integrated coherent
  cross section for Z 92–98 is 1–34% lower, most just *below* K edges — the anomalous-scattering pattern —
  so the integrated cross section may already carry it and only the angular sampling (F² alone) lack it.
- **Shows up when:** sub-30 keV high-Z work; accurate Rayleigh-peak shape.

---

## Validation bookkeeping

### GEANT4 references are still scored at the old 1.5 keV FEP window `[medium — validation]`

**Effect.** `physics/FepWindow.h` consolidated "what counts as full-energy" onto one constant,
`kDefaultFepWindowKeV = 0.75` keV. Every committed GEANT4 reference in
`tests/data/geant4_reference/` predates that and was generated at a 1.5 keV half-window, as were the
several hundred numeric expectations in the C++ suite. Anything compared against them therefore pins
1.5 keV instead of the default: `kTestFepWindowKeV` (`tests/test_fep_window.h`, used by the suite and
`cascade_ref_common.h`) and `kReferenceFepWindowKeV` (`examples/benchmark_mc_configs.cpp`, the
writer of the `our_*_multi.csv` input to the `profiling/compare_validation.py` gate). Until the pins go, the gate and the tests are not
measuring the library's own default. Two producers deliberately follow the new default instead and
so are NOT comparable to the committed references: the harness itself
(`EventAction`, defaulted, no CLI override) and `tools/geant4_validation/generate_all_spectra.cpp`. The gap is real where near-peak forward scatter is: config 8 at
59.5 keV gives FEP 7.479e-02 +/- 0.161% at 1.5 keV vs 7.248e-02 +/- 0.163% at 0.75 keV, -3.09% at
z = -13.5. Bare-NaI config 1 at 60/662/1332 keV is unaffected (|z| <= 1.0 at 0.1% precision).

**Implementation.** Settle the default first (0.75 keV suits HPGe; the MC has no resolution
smearing, so the window is really a choice about how much near-peak scatter to credit). Then
regenerate the GEANT4 references with `EventAction::SetFepWindowKeV` set to the same value — the
harness setter exists but has no CLI flag yet, so add one to `main.cc` — re-baseline the suite's
expectations, and delete both pins plus the note at the top of
`tests/data/geant4_reference/README.md`.

**Shows up.** Any comparison of CeeLo-at-default against the committed references; low-energy and
scattering-heavy configs (8, 7) most of all.

## Error estimation

### D3 — no systematic-error budget `[low]`
- **Effect:** reported σ is purely statistical; documented model biases (LaBr3 ≥1 MeV residual, cfg-8
  59 keV near-peak scatter, cfg-12 59 keV data difference, skin-escape; the former ≥2 MeV FEP family and
  CZT escape were fixed Aug 2026) exceed σ_stat in their regimes, so a user computing
  z = Δε/√(σ_MC²+σ_G4²) gets spuriously large |z| that is a known model effect, not a discrepancy to chase.
- **Implementation:** document (and optionally expose) an energy/config-dependent σ_sys so z-scores in the
  bias-dominated regimes use √(σ_stat²+σ_sys²); at minimum annotate the affected energies.
- **Shows up when:** computing MC-vs-G4 z-scores in the known-bias regimes.

---

## Source geometry (new shapes)

Current shapes: point, solid/annular cylinder, rectangular box, solid/hollow sphere, Marinelli beaker
(`SourceGeometry::Shape`). Each supports a source material/density (self-attenuation), uniform (or
exponential-depth) activity, isotropic emission, arbitrary rotation/off-axis placement, and full-MC
source transport. Sphere/cylinder/box ray-intersection primitives exist (used for shield shells). A
thin disk / planar deposit is **already expressible** as a degenerate thin solid cylinder
(uniform-volume sampling and ray-cylinder self-attenuation are both correct in that limit) — only a
convenience constructor and an optional "massless deposit" (self-attenuation off) flag are missing,
not a new shape.

(Done — **solid + hollow sphere** (`set_spherical_source`) and **hollow/annular cylinder**
(`set_cylindrical_source(..., inner_radius)`) landed with the cascade geometry-parity work; see
DESIGN.md "Source shapes" and `tests/test_spherical_source.cpp`. **Hollow rectangular box shell**
(`set_rectangular_source(..., inner_half_dims)`, crate/container walls) landed with the Stage-B
validation-matrix work; see `tests/test_hollow_box_source.cpp`.)

### Position bias ignores hollow-source voids `[low — coverage]`
- **Effect:** `sample_biased_position` / `compute_auto_bias_params` ignore `inner_radius` /
  `inner_half_dims` (they sample the full solid proposal with solid-volume weights), so enabling
  position bias on a hollow source emits from the void with wrong weights. Unbiased sampling
  (the default for these shapes) is correct.
- **Implementation:** reject-and-reweight (constant factor V_solid/V_shell) or annulus-aware lateral
  sampling in `sample_biased_position`, or assert position bias off for hollow sources.
- **Shows up when:** hollow sources large enough that auto position bias would otherwise engage.

(Done — the "pipe wall source with attenuating contents" half of this item is resolved by
`add_source_core()`: a hollow source's interior can now be filled with attenuating layers instead of
being a mandatory void. See DESIGN.md "Source shapes" and `tests/test_source_core.cpp`.)

### Extended-source TOTAL has sub-percent structure at the ends of the range `[low — fidelity]`
- **Effect:** with source-electron transport enabled and both sides at high precision, the
  uncored soil-shell control (G-C) total sits **+0.62% (z 2.1) at 238.6 keV, +0.74% (z 2.9) at
  583.2, −0.05% at 911.2 and −0.63% (z −3.6) at 2614.5**. Sub-1% everywhere and inside the ~1%
  goal, but now statistically resolved rather than lost in noise. The cored geometries (G-D/G-F/
  G-G) show the same size and sign, so it is a property of the extended-source path, not of cores.
- **Implementation:** unknown; candidates are the secondary-channel normalization above the pair
  threshold (which now slightly over-delivers at 238-583 keV while still under-delivering at
  2614) and the scattered-stream placement documented for config 8.
- **Shows up when:** comparing extended-source TOTAL efficiency at high statistics.

### GDML export emits no rotation for source volumes `[medium — silent mismatch]`
- **Effect:** `write_gdml()` places every source volume with a `<position>` and **never** a
  `<rotation>`, so a `CylinderSideOn` source (90° about y) or any rotated box is exported
  **unrotated**. The MC honours the rotation, so the exported scene is not the simulated scene —
  silently, and only for configurations no current benchmark uses.
- **Implementation:** emit a `<rotation>` (or `<rotationref>`) alongside the position for the
  outermost source volume, taking the angles from `cyl_rotation_` / `rect_rotation_`. Nested
  daughters share the mother's frame, so only the world placement needs it.
- **Shows up when:** validating a side-on cylinder or a tilted box source against GEANT4.

### add_core() overflow is silently clamped `[low — API]`
- **Effect:** core thicknesses summing to more than the source's inner extent are clamped to zero
  extent (`.cwiseMax(0.0)` in `rebuild_layer_dims`) rather than rejected, so an over-specified stack
  quietly loses its innermost layers. The hollow-source precondition IS asserted; the total is not.
- **Implementation:** assert (or throw) when the accumulated core thickness exceeds the cavity.

### Source-core validation gaps `[low — validation coverage]`
- **Effect:** three things the Sep 2026 source-core campaign did not pin.
  (a) A **partly-filled core** — cores that do not reach the centre, leaving a real non-attenuating
  cavity — has never been GEANT4-validated; it is covered only by unit tests, although the design
  doc called it out explicitly. (b) The **GPS-confinement ratio control (G-E) exists only for the
  sphere**; the cylinder and box rely on their direct CeeLo-vs-G4 agreement instead. (c) The
  **partly-filled and FAR (~50 cm) cored cases have no G4 anchor at all**, so the far-field
  cone/two-stream path for a cored source is unvalidated.
- **Implementation:** a G-D-like case with a half-filled cavity at NEAR, and G-E-style ratio
  controls for the cylinder and box, would close (a) and (b) for a handful of 16M G4 runs.
- **Shows up when:** a source whose cores leave a central void, or a cored source at ≫10 cm.

### Extended-source TOTAL is a coherent family, not just noise `[low — fidelity]`
- **Effect:** across the 16 cored + control rows, FEP is statistically clean (χ²/n ≈ 1.1) but
  **TOTAL is not: χ²/n ≈ 3.2**, with 5 rows at |z| > 2 against ~0.7 expected and a positive mean.
  Cored and uncored geometries show the same sign and size, so it is a property of the
  extended-source path rather than of cores. Every deviation is still sub-1%.
- **Implementation:** unknown; see the neighbouring item on sub-percent structure at the ends of
  the range, which this generalises.
- **Shows up when:** comparing extended-source TOTAL efficiency at high statistics.

### Semi-infinite half-space (wall / floor / soil slab) `[medium — coverage + estimator]`
- **Effect:** large-area surface or volumetric sources (contaminated wall/floor, bulk soil/concrete).
  Approximable today as a **large, deep** rectangular box (thin box only for surface contamination;
  volumetric needs depth of several mean-free-paths), but uniform-volume box sampling wastes nearly all
  histories far off-axis — correct but impractically slow without `lambda_lateral` position biasing.
- **Implementation:** a first-class half-space template with lateral importance weighting baked in
  (instead of hand-tuned `lambda_lateral`), analytic lateral truncation (no "how big is big enough"
  guesswork), per-area (Bq/cm²) / per-mass (Bq/g) normalization, and possibly an angular-flux estimator
  rather than brute-force volume sampling. Air-path effect (below) matters more here than for compact
  near-contact sources.
- **Shows up when:** in-situ / large-area environmental measurement with standoff.

### Lower-priority shapes (need new intersection / sampling code) `[low — coverage]`
- **Effect:** truncated cone / frustum (funnels, tapered vials, conical piles); bottle / cylinder-with-
  neck (partial fill + headspace, where the fill line and neck matter); multiple spatially-separated
  active regions in one run; non-uniform activity beyond exponential-depth (radial gradients, surface
  plate-out, shell-weighted). None are expressible with the current four shapes.
- **Implementation:** frustum needs a new cone-intersection routine; bottle is a composite solid; multi-
  region needs the source loop to sum over regions; non-uniform activity extends the sampler weighting
  (the `PositionBiasConfig` / `DepthDistribution` hooks already exist) rather than the geometry.
- **Shows up when:** specialized sample containers; composite / graded sources.

---

## Detector geometry

### Germanium shows a ~0.1% energy-structured offset vs GEANT4 at 60–122 keV `[low — residual]`
- **Effect:** on the HPGe pair, CeeLo sits +0.135% above GEANT4 in FEP at 88 keV (z ≈ 3.5–3.9 at
  ~0.03% MC precision) and about −0.1% at 122 keV, with ~+0.05% at 59.5 and ~−0.03% at 45. The
  offset is **the same for the sharp and bulletized crystals** (88 keV: cfg 25 +0.135%, cfg 26
  +0.136%, difference +0.000 pp), so it is a property of germanium transport, not of the fillet —
  it cancels in the bulletization comparison. Originally mistaken for a bulletization residual;
  higher CeeLo statistics showed the 45 keV part was scatter and the rest was common to both.
- **Implementation:** the shape — small, energy-structured, sub-percent — matches the EPICS2023 vs
  GEANT4-EPICS2017 photon-evaluation difference already established for Fe/Cr/Ni (commit 1ad196c).
  Check it the same way: compare µ/ρ(Ge) from `element_data.cpp` against NIST XCOM (Hubbell &
  Seltzer, SRD 126) over 40–150 keV, and against G4's own `EmCalculator` export via
  `g4_extract_xs`. If EPICS2023 is closer to NIST, this is CeeLo being more accurate and the item
  becomes a documented evaluation difference rather than an accuracy gap.
- **Shows up when:** HPGe below ~150 keV. Well inside the ~1% goal.

### Dead layer is not exported to GDML `[low — validation coverage]`
- **Effect:** `write_gdml()` writes the full crystal as `active_crystal` and only *notes* the dead layer
  in a comment, so any GEANT4 comparison of a dead-layer geometry is really comparing different solids.
  Configs 25/26 therefore carry no dead layer, and the bulletized *active* fillet (the inward-offset arc
  in `trace_cylinder_geometry`) has unit-test coverage but no G4 anchor.
- **Implementation:** export the crystal as two nested solids and score only the inner one — the outer
  shell is the same material, so it must be a separate logical volume rather than a subtraction. The
  harness already keys on the `active_crystal` volume name, so only the export and the shell's naming
  change.
- **Shows up when:** validating HPGe with realistic 0.5–1 mm dead layers, especially below ~100 keV.

### ANGLE `.outx` detector definitions are not importable `[low — convenience]`
- **Effect:** the crystal geometry CeeLo now models (radius, length, `bulletizingRadius`, bore radius /
  depth / `rounded`, inactive-Ge thicknesses, endcap and housing layers) is exactly what an ANGLE `.outx`
  file already carries, but it has to be transcribed into `set_*` calls by hand.
- **Implementation:** a small reader mapping `<crystal>`/`<core>`/`<inactiveGe>`/`<endCap>`/`<housing>`
  onto `set_detector` + `set_bullet_radius` + `set_bore_hole` + `set_dead_layer` + `add_attenuator`.
  Note `.outx` is in mm. Most naturally lives on the InterSpec side rather than in the library.
- **Shows up when:** reproducing a vendor-characterised HPGe.

---

## Other known limitations (documented; tracked for awareness)

### High-coverage partial graphs can retain small topology-free residuals `[low — nuclear-data policy]`
- **Effect:** coverage is measured in repaired transition-occurrence weight, including memberless E0
  and highly converted transitions. Every accepted graph has at most 0.01 absolute residual occurrence;
  below 99% relative coverage it additionally requires complete RDM starts that replay the matched
  transitions. The remaining small topology-free residuals are retained as independent categorical
  transitions because exact correlation information is absent.
- **Known cross-evaluation case:** high-coverage Th-234 remains on the inferred DAG. Its old RDM exact
  starts disagree with the newer evaluated 63.29/92.8-keV photon yields; forcing those starts would move
  the result away from the photon evaluation the cascade currently models.
- **Implementation:** replace the independence approximation only if a future enrichment maps the
  residual transition topology, or supplies authoritative joint correlations for it.
- **Shows up when:** unmatched transitions carry less than one percent of the branch occurrence sum.

### Air gap not simulated `[documented]`
- **Effect:** air between the outermost source shield (or the source) and the detector face is not in the
  MC transport; intended to be applied later as a deterministic analytic correction.
- **Implementation (investigate):** split the air correction by component so it stays cheap. For the
  **unscattered / full-energy** component, apply an analytic (analog) transmission T_air = exp(−μ_air·ℓ)
  along the source→detector path — exact, essentially free, removes the build-out of FEP at low energy /
  long standoff. For the **scattered** component (total efficiency / spectrum: small-angle air scatter
  that still reaches the crystal, plus scatter-out), account for it with a GadrasShieldScatter-style
  analytic scatter kernel rather than tracking air explicitly in the MC — or some other deterministic
  add-on — so the bulk MC speed is preserved. Validate the split against a G4 run that includes the air
  volume.
- **Shows up when:** long source-detector standoff, low energy; matters most for the half-space / in-situ
  geometry above.

### Doppler-broadening approximations `[documented; NOT an FEP/total effect]`
- **Effect:** the F(p_z) profile-reweighting factor is omitted (few-% tail) and there is no relaxation
  cascade from the Compton vacancy. Compton-edge **shape** only — never cite for FEP/total discrepancies.
- **Shows up when:** detailed Compton-edge shape comparisons.
