# CeeLo — Open Items / TODO

The single register of known shortcomings, deferred features and fidelity/speed work.
**When an item is resolved, delete it in the same commit that fixes it** (and update
DESIGN.md if it described the gap).

Each item gives its **Effect** (what is wrong or what it buys), **Implementation** (what it
would take) and **Shows up** (the regime where it matters), plus a severity tag relative to
the ~1%-vs-GEANT4 goal. Areas are listed in priority order and items within an area are
sorted the same way; "Next up" is the short cross-area ranking. Items measured to be
sub-statistical in the gated configs are kept for awareness.

## Next up

1. [Wire the GEANT4 gate into ctest](#compare_validationpy-is-not-wired-into-ctest) — cheap; stale
   CeeLo CSVs have hidden real drifts and invented phantom ones.
2. [Load-robust precision tests](#precision-targeted-mc-tests-are-wall-clock-capped) and the
   [IC-electron sum-peak test](#ic-electron-sum-peak-test-asserts-an-unresolvable-ordering) — flaky
   tests erode trust in the suite.
3. [Regenerate the GEANT4 references at the 0.75 keV window](#geant4-references-are-still-scored-at-the-old-15-kev-fep-window)
   — until then the gate does not measure the library's own default.
4. [Rejected-branch pair links](#rejected-branch-pair-links-cannot-all-be-honored) and
   [infeasible ICC/intensity pairs](#infeasible-iccintensity-pairs-can-cost-branches-their-level-scheme)
   — summing for U-235, Pu-241, Sb-125 silently switches model.
5. [GDML export: source rotation](#gdml-export-emits-no-rotation-for-source-volumes) and
   [dead layer](#dead-layer-is-not-exported-to-gdml) — silent validation mismatches.
6. [Air between source and detector](#air-between-source-and-detector) — validate the analytic
   attenuation and decide its default.

---

## Validation and test infrastructure

### `compare_validation.py` is not wired into ctest
`[medium — validation]`
- **Effect:** the gate compares the GEANT4 references against whatever `our_*_multi.csv`
  sits in `build/examples/` (written by `benchmark_mc_configs --config N --precision p`;
  `--mc-dir` overrides), and nothing checks that those files come from the current engine.
  When the engine moves and they are not regenerated, the gate reports staleness as physics
  (a "low-energy Marinelli deficit" and a "coherent ≥2 MeV FEP family" were both stale
  files), or misses real drift. The cascade block has the same problem
  (`cascade_observables` → `our_cascade_multi.csv`). Separately, not every reference in
  `tests/data/geant4_reference/` records its GEANT4 physics list in its header.
- **Implementation:** a ctest case that regenerates the CeeLo side at modest precision and
  compares, or at least a staleness check (the CSVs record their producing commit and the
  test fails when `src/` is newer); include the cascade block. New references stamp their
  physics list.
- **Shows up:** after any engine change that lands without regenerating the CeeLo side.

### Precision-targeted MC tests are wall-clock-capped
`[medium — test design]`
- **Effect:** ~15 test files set `target_fep_rel_precision` together with
  `termination.max_wall_seconds` (30–120 s), and whichever comes first wins. On a loaded
  machine a run can stop on wall time with the target unmet while its assertions use fixed
  tolerances, so a truncated run fails on statistics with no sign that is what happened.
  The tests seen flaking under load (`test_asymmetric_shield`, `test_position_bias`,
  `test_marinelli_beaker`'s `source_electron_transport_increases_total`) are all in this
  set; the last compares two independent 2%-precision runs against a 5% bound (~1.8σ even
  at target). The mechanism is the leading suspect, not a demonstrated cause.
- **Implementation:** first, on failure report `StopReason` and the achieved
  `fep_uncertainty` (both returned by `compute()`), so the next failure arrives as evidence.
  If truncation is confirmed, assert on a z-score from the returned uncertainty instead of a
  fixed tolerance across the `precision_config` family.
- **Shows up:** ctest on a busy machine, especially beside a build.

### IC-electron sum-peak test asserts an unresolvable ordering
`[low — test design]`
- **Effect:** `test_cascade_summing::ic_electron_g4_magnitude_ba133_2cm` requires the
  Ba-133 `sum356+81` coincidence peak to land closer to GEANT4 with IC electrons ON than
  OFF. ON/G4 = 1.009 and OFF/G4 = 0.998 against a reference with its own ~1% statistics,
  so which is closer is a coin flip — the test can fail at random. The resolvable quantity
  is the sign of the ON-vs-OFF shift (a sum peak gains counts when singles lose them).
  `feature/FullSpecApprox` c6b1436 already makes that change.
- **Implementation:** port c6b1436 (assert the shift's direction, bound the magnitude by
  the ON/G4 band). Making the ordering itself testable needs a ~25× larger Ba-133 GEANT4
  ion run.
- **Shows up:** intermittent `test_cascade_summing` failures.

### GEANT4 references are still scored at the old 1.5 keV FEP window
`[medium — validation]`
- **Effect:** `kDefaultFepWindowKeV` is 0.75 keV, but the GEANT4 references for configs 1–12
  and 25/26 were generated at 1.5 keV, as were the C++ suite's numeric expectations, so the
  gate and tests pin 1.5 keV (`kTestFepWindowKeV` in `tests/test_fep_window.h`, used by the
  suite and `cascade_ref_common.h`; `kReferenceFepWindowKeV` in
  `examples/benchmark_mc_configs.cpp`). Until the pins go, they do not measure the
  library's default. Configs 27/28 are already at 0.75 keV on both sides; the harness and
  `generate_all_spectra` score at 0.75 keV, so new harness runs are not comparable with the
  old references. The gap is real where near-peak forward scatter is: config 8 at 59.5 keV
  differs by −3.09% (z −13.5); bare NaI does not move.
- **Implementation:** settle the default (the MC has no resolution smearing, so the window
  is a choice of how much near-peak scatter to credit); add a `--fep-window` flag to the
  harness `main.cc` (FS has one), regenerate the references, re-baseline the suite, and
  delete both pins and the note in `tests/data/geant4_reference/README.md`.
- **Shows up:** comparing CeeLo-at-default against the references; configs 7 and 8 most.

### Configs 20/21 (Pb/W-wall boxes) have no GEANT4 reference on main
`[low — validation coverage]`
- **Effect:** the high-Z wall boxes are exportable but ungated here; a suspected +2–3.7%
  total offset at 2614 keV was never anchored. `feature/FullSpecApprox` has 16M-event
  option4 references at 662/2614 keV and gate rows for both: they show no common offset
  (Pb −0.1%, W +1.45% at 2614 keV before the crystal-walk fixes; ≤ 1.7% now), with W's
  662 keV FEP the worst row.
- **Implementation:** port the two reference CSVs and their `compare_validation.py` rows
  (with that branch's tightened cfg 5 tolerance) to main.
- **Shows up:** Z ≳ 74 wall boxes.

### Fully independent Conditional summing-factor gate
`[low — validation]`
- **Effect:** the cascade gate (`tests/data/geant4_reference/cascade_summing_multi.csv`,
  shared by `compare_validation.py` and ctest) checks Conditional rows only as
  `A_FR·k_cond/k_full`, i.e. anchored to FullRealization's per-decay normalization; a truly
  independent check needs a GEANT4 "no-summing" baseline.
- **Implementation:** time-stagger in the harness — in a
  `G4UserTrackingAction::PreUserTrackingAction`, offset each genuine decay gamma (creator
  process `RadioactiveDecay` only; offsetting every gamma would split a primary's own
  fluorescence/brems from its FEP) by `TrackID × 1 ms`, then bin `(t_prestep, edep)` into
  ms windows: one entry per window is the no-summing spectrum, their sum the with-summing
  one, so `k_G4 = A_corr / A_nosum` from a single run (`--no-summing`). β⁺ pairs and
  coincident x-rays need an extra lineage rule; Co-57/Co-60/Ba-133 are clean.
- **Shows up:** wanting an estimator-independent Conditional check.

---

## Monte-Carlo sampling (cone)

Found by review on 2026-11-03 and measured there with direct Monte Carlo; not yet re-verified.

### Near-axis cone shortcut samples about ±z
`[medium — accuracy]`
- **Effect:** `sample_cone_about` (`EfficiencyCalculator.cpp`) samples about ±z whenever the
  wanted axis is within 0.81° of it (`cos_rot > 0.9999`), so the cone is tilted off the detector.
  For a 3x3 NaI with the source at (0.5, 0, −40) cm, the coned total efficiency was 1.7292e-3 ±
  7.8e-7 against 1.7585e-3 ± 6.6e-6 isotropic: −1.67%, z ≈ 4.4.
- **Implementation:** an exact rotation. Keep the identity branch for a source exactly on the
  axis, so those runs stay bit-identical.
- **Shows up:** point sources slightly off the axis, far away.

### The legacy cone half-angle misses thin detectors
`[medium — accuracy]`
- **Effect:** `compute_raw_cone_half_angle`'s six candidate points miss the widest angle for thin
  crystals (L/R about 0.5).
  - Example: a planar Ge (R 2.55, L 1.25, 0.05 Al) with the source at d 5.2 cm, θ 86°. Coned vs
    isotropic: FEP −2.1% (z ≈ 4.2), total −2.7% (z ≈ 13.5).
  - Up to 10% of the directions that hit the detector fall outside the padded cone at d = 2a.
  - The side & back stage uses this cone for sources beyond z_max. For thin detectors that costs
    about 0.3% at d_far and up to 1.1% close behind.
- **Implementation:** use the exact `envelope_bounding_cone` wherever bit-identity with old runs is
  not required, e.g. whenever general cone sampling is on.

### `Geometry::outer_z_extent()` does not accumulate front thicknesses
`[low — accuracy]`
- **Effect:** with stacked layers (all `z_start` = 0; a 0.4 cm gap front, then a 0.1 cm can
  front), it reports z_min = −0.4. The traced front surface and `outer_envelope()` are at −0.5.
  Both the legacy and the general cone use the smaller extent (`RayTrace.cpp`).
- **Implementation:** accumulate as `outer_envelope()` does.

### General cone sampling re-derives the cone every history
`[low — performance]`
- **Effect:** with `enable_general_cone_sampling` on, `sampling_cone` calls
  `envelope_bounding_cone` (about 19 µs CPU) once per history. That happens on the per-position
  cone path and on the source-effects direct-stream and mixture paths.
  - For an extended source beside the detector, this can cancel the roughly tenfold saving in
    histories.
  - The side & back generation stage uses point sources, so it is unaffected.
- **Implementation:** cache the cone per source position bin, or bound it once over the source
  volume.

---

## Photon transport physics

### Air between source and detector
`[medium — coverage]`
- **Effect:** air between the outermost source layer and the detector is not transported.
  `set_air_attenuation(AirAttenuation::AnalyticNoScatter)` (opt-in, default `None`)
  multiplies each event's weight by exp(−μ·d) over the source-to-detector air path (μ
  without Rayleigh), which covers the unscattered peak, but it has no test and no GEANT4
  comparison; scatter-in and scatter-out are not modelled.
- **Implementation:** test and GEANT4-validate `AnalyticNoScatter` (an air world, long
  standoff, low energy) and decide whether it should be the default; handle the scattered
  component with a deterministic kernel rather than tracking air in the MC.
- **Shows up:** long standoff at low energy; half-space / in-situ geometries.

### cfg 8 Marinelli FEP −1.6% at 59 keV
`[low — residual]`
- **Effect:** the water Marinelli runs FEP −1.6% (z −5.5) at 59 keV and −0.8% at
  100–662 keV; totals ≤0.6%. The GEANT4-harness u/s decomposition localises it: the
  unscattered FEP `eps_u` matches to 0.01–0.3% at every energy, while the small scattered
  FEP stream `eps_s` (forward Rayleigh / small-angle Compton landing in the window) is ~15%
  low. Behind deep iron (configs 27/28) the scattered share of the 60 keV peak agrees to ~2%
  (30.0% vs 30.6%), so the deficit does not generalise to Rayleigh-dominated cases.
- **Implementation:** improve near-peak forward-scatter recovery at low energy, possibly
  with a crystal-escape albedo-recovery model (GEANT4 also recovers crystal-escaped e⁻/γ
  backscattering off the surrounding water; CeeLo's re-entry is ≤3-bounce photon-only).
- **Shows up:** Marinelli / wrap-around extended sources near 60 keV (FEP only).

### LaBr3 high-energy FEP residual
`[low — residual]`
- **Effect:** config 3 (LaBr3 2"×2") FEP sits +0.72% (z 2.3) at 2000 and +0.63% (z 1.9) at
  3000 keV in the gate (+0.63% z 4.0 / +0.49% z 2.6 in a dedicated 0.16%-precision study);
  NaI configs read up to +0.5% at 2–3 MeV.
- **Implementation:** candidates are Bohr straggling variance slightly large for LaBr3 (no
  shell/density corrections) or a cfg-3 reference offset.
- **Shows up:** LaBr3 photopeaks ≥2 MeV, at the +0.5% level.

### Pb-attenuator total-efficiency excess
`[low — residual]`
- **Effect:** config 7 (2 mm Pb attenuator) total +0.21 ± 0.08% / +0.19 ± 0.07%
  (z +2.8 / +2.5) at 2000/3000 keV; FEP agrees. It read +0.25 / +0.27% before a positron
  that leaves the geometry stopped making a 511 keV pair. On feature/FullSpecApprox, whose
  electron walk differs, the same fix took it to +0.09% / +0.06%.
- **Implementation:** untriaged; plausibly the attenuator's pair-production / secondary
  treatment.
- **Shows up:** high-Z attenuators ≥2 MeV, totals only.

### EADL K line energies are off by up to 0.3 keV for heavy elements
`[low — spectrum fidelity]`
- **Effect:** the K line energies in `relaxation_epics_data.cpp` are EADL subshell
  binding-energy differences. For Ge, Cd, Te and I they sit within 0.04 keV of the
  measured lines, but Pb Kα1 is 75.25 keV against 74.97 keV measured (Kα2 73.04 against
  72.80); GEANT4 option4 puts Pb Kα1 at 74.97. Now that source and shield K x-rays are
  emitted, a Pb shield's x-ray peaks on HPGe sit 0.28 keV high (the FWHM there is ~1.1 keV).
- **Implementation:** replace the line energies with measured transition energies
  (Deslattes et al., *Rev. Mod. Phys.* **75** (2003) 35; NIST SRD 128) in
  `generate_eadl_relaxation.py`, keeping EADL's intensities and yields; regenerate.
- **Shows up:** Pb, W, U and Th x-ray peaks on HPGe.

### Pair production near threshold is unarbitrated
`[low — residual]`
- **Effect:** EPICS2023's pair cross section is below GEANT4's parameterization near
  threshold: at 1332 keV Ge −9%, Cd −5%, I −3%, Na −6%, Al −6%; at 1.2 MeV I −21%; at
  2614 keV within ±1% except Pb +2.3%. In GEANT4 pulse-height comparisons the single and
  double escape peaks at 1332 keV read 5–13% low accordingly (NaI SE −5.4%, HPGe SE −13%,
  CZT DE −7%), and the 511 keV line behind Pb shells at 2614 keV a few percent high.
- **Implementation:** compare both against NIST XCOM pair production (nuclear + electron
  field) at 1.1–3 MeV for Na, Al, Ge, Cd, I and Pb, as was done for the iron
  photoelectric difference; document it in DESIGN.md as an evaluation difference, or
  correct the data if EPICS2023 is the outlier.
- **Shows up:** single/double escape peaks below ~2 MeV; pair-fed 511 keV lines.

### L/M-shell fluorescence x-ray escape not modelled
`[low]`
- **Effect:** only K-shell fluorescence is transported; L x-rays (e.g. Pb Lα ≈ 10.5 keV) are
  deposited locally, over-predicting local deposit and missing L-escape features. On a bare
  1 cm CZT, GEANT4 shows Cd/Te L-escape lines at E₀ − 3.1…4.3 keV (~0.02% of the FEP) that
  CeeLo lacks, and a bare self-attenuating thorium sphere reads ~1% low in total, plausibly
  its L x-rays.
- **Implementation:** emit L lines from the EADL data already in
  `relaxation_epics_data.cpp`.
- **Shows up:** thin high-Z detectors (CZT), PE in a high-Z attenuator skin, low-energy
  spectra.

### Pair-production energy sharing and transport above 3 MeV
`[low — fidelity, coverage]`
- **Effect:** the crystal splits the pair kinetic energy exactly 50/50
  (`PhotonTransport.cpp`, `ke_each = 0.5 * deposited_local`) instead of sampling the
  Bethe–Heitler sharing, so both leptons ride the same deterministic CSDA step grid; it
  leaves a 46.2 keV comb in pair-history deposit spectra at 4–9 MeV (seen in full-spectrum
  truth) and is a candidate for small escape-peak shape errors. The source side already
  samples the split. More broadly, nothing above 3 MeV is GEANT4-gated although the photon
  tables reach 20 MeV, and photonuclear absorption is absent (as in EPDL/XCOM; a few
  percent of Pb/U attenuation at 12–15 MeV).
- **Implementation:** sample the split from the screened Bethe–Heitler cross section and
  A/B it; a GEANT4 gate row at 10–15 MeV for a bare and a Pb-shielded detector;
  photonuclear only if it ever matters.
- **Shows up:** escape-peak/continuum shape above ~3 MeV; reaction gammas (e.g.
  N-14(n,γ) 10.8 MeV).

### cfg 7 at 100 keV (Pb K edge)
`[GEANT4-side; excluded from the gate]`
- **Effect:** ~−25 to −31% total at 100 keV, traced to GEANT4's in-tracking μ_total(Pb)
  being smeared across the K edge (DESIGN.md). Separately, CeeLo's escaped Pb Kα/Kβ ratio is
  0.67 against GEANT4's 0.82 (a simple emission × attenuation estimate agrees with GEANT4).
- **Implementation:** none on CeeLo's side for the edge; verify CeeLo's Pb K-line branching
  ratios and energies against EADL.
- **Shows up:** just above high-Z K edges behind thick high-Z attenuators.

---

## Cascade / true-coincidence summing

### Rejected-branch pair links cannot all be honored
`[high — decay data]`
- **Effect:** branches whose level scheme is rejected fall back to flat member marginals
  plus directional pairwise links, which are not jointly consistent (a whole-database census:
  10,355 linked pairs whose implied joint exceeds `min(P(a),P(b))`; 3,172 of 131,851 pairs
  whose two directional reconstructions disagree by >1%). FullRealization projects them onto
  a deterministic one-parent Bayesian forest (every pair joint Fréchet-bounded, every
  marginal preserved), so incompatible links go unused; Conditional and analytic do not
  invent gamma-gamma sum-fed joints for rejected branches. All three set per-peak
  `summing_model_complete=false` when an invalid primary or a possible invalid-branch photon
  pair affects the peak.
- **Implementation:** recover authoritative level topology or a consistent multivariate
  decay law for the rejected records. Averaging directional conditionals or clamping pairs
  independently cannot define a coherent 3+-member realization. Analytic/Conditional should
  enumerate the forest's pair joints only if that projection becomes part of their contract.
- **Shows up:** any peak marked `summing_model_complete=false`, especially U-235, Pu-241,
  Sb-125.

### Infeasible ICC/intensity pairs can cost branches their level scheme
`[medium — decay data]`
- **Effect:** some records pair a gamma intensity with an ICC for which
  `I_gamma*(1+alpha) > 1`. The adapter bounds the transition and rejects a graph whose
  reconciled feeding stays infeasible; rejected branches use the pairwise model above. The
  regenerated-XML audit reports 4566 valid, 754 raw-only, 1782 partial-only and 1328
  raw+partial branch instances, with 323 E0 repairs and 804 intensity caps in accepted
  graphs (a distinct parent→child recount is still needed before quoting a rejection
  rate). Rejected examples with T½ > 1 y: U-235 → Th-231 (raw_total_feed 2.37), Th-229 →
  Ra-225, Hf-172 → Lu-172m, Pu-241 → U-237, Pa-231 → Ac-227, Sb-125 → Te-125, plus Am-242m,
  Cf-249, Bk-247, Cm-243, Np-235, Ho-166m, Po-208, Sn-126, Es-252; also Pa-234 → U-234
  (1.25). U-235's summing factors change 12–850× across the switch, so no absolute value
  from a rejected nuclide's level path is trustworthy. The GEANT4-anchored cascade gate
  (Am-241, Ba-133, Co-57, Co-60, Eu-152, I-125, Na-22) cannot see a regression here.
- **Not a matching failure**: U-235's 19.59 keV gamma is correctly matched (Th-231 level
  5 → 4, alpha = 114.7); the input evaluations disagreed, and the regenerated XML carries the
  curated ENSDF override (0.00583 per decay).
- **Already done:** a member's gamma intensity is capped at `1/(1+alpha)`, the ceiling that
  holds whichever number is wrong (I-125's Te-125 transition loses 0.7%; U-235's 19.59 keV
  loses a factor of 70, removing false sum peaks at 205.31 and 221.40 keV).
- **Implementation:** uncurated records need a per-transition ENSDF check of whether the
  tabulated intensity is the gamma or the transition rate (435 of the historical 1052
  violations sit within 2% of the bound, where the two readings differ 15×); no heuristic
  separates them.
- **Shows up:** U/Pu-chain nuclides with strongly converted low-energy transitions; the
  calibration set is unaffected.

### Analytic cascade summing: remaining approximations
`[low]`
- **Effect:** `compute_cascade_analytic` matches FullRealization to ≲1% at all geometries
  including contact for the mainstream γ-γ / EC-x-ray / β⁺ cases (Co-57, Ba-133, Na-22,
  Am-241, Eu-152, Co-60; |z| ≤ 0.97 at 2 cm). Known residuals:
  - **Escape-peak recovery** (`Π(1−p·ε_tot)` treats the peak as a delta): misses events where
    the primary's escape peak plus a coincident x-ray sum back into the window — I-125 35.49
    @2 cm −4%, Ba-133 81 @2 cm −6%, Eu-152 121.78 −1.5%; ~0 by 20 cm. Not fixable from
    ε_FEP/ε_tot alone; these rows are report-only in the gate.
  - **W(θ)** uses the collinear limit `g = W(0)` for correlated summing-in pairs (recovers
    Co-57 122+14 → 136 to ~0.2%); slight over-correction possible at very close geometry.
  - **Triple-fed summing-in** is enumerated only for gamma lines on distinct transitions;
    x-ray / EC / annihilation triples, 4+-member sums and topology-free residual gammas are
    omitted (ε³-small, and a perf necessity).
  - **Non-gamma primaries**: a 511 primary uses a marginal-product fallback; a peak at a
    pure vacancy x-ray energy returns `found=false`.
  - IC-electron deposition is whatever the ε provider models; Kα→L secondary vacancies are
    dropped (summing-out is idempotent, summing-in ε²-small).
- **Implementation:** a deposit-spectrum-aware provider (escape energy + fraction); a real
  angular joint acceptance; EC/annihilation triples if a case needs them.
- **Shows up:** very low-energy lines just above the detector K edge at contact;
  sub-0.5% work on x-ray-rich EC nuclides at contact.

### Transition-group (flat) x-ray model: remaining approximations
`[low]`
- **Effect:** the flat model (`vacancy_xray_model=false`; InterSpec's reference lines and
  sum-peak check) carries K–K pairs (DESIGN.md) but: x-rays and gammas are linked as
  independent (a Pb-203 K–K pair needs the 279 keV transition converted, yet the link treats
  that gamma as coincident; K–K sum survival reads ~10–20% low at contact); EC capture and
  conversions are independent in `D2` (+0.2% Pb-203, +2.8% Ba-133, +5.9% Te-123m); the MC
  ignores `self_pair_prob`; no L x-ray pairs or K→L relaxation; and gamma-gamma links are
  SandiaDecay's coincidence fractions verbatim, which disagree with the level scheme for
  ~20k of ~27.5k links (Lu-177 all ~0.1×; Eu-152 125.7 → 562.9 keV 0.37 vs 0.99; clean:
  Ba-133, Co-60, Ra-226, Tl-208, I-131).
- **Implementation:** the vacancy model is exact but too slow for InterSpec's redraw
  (36–78 ms per window for Eu-152); cheaper: condition the flat x-ray→gamma links on the
  vacancy source. The gamma-gamma links are a decay-data question first.
- **Shows up:** x-ray–gamma or x-ray–x-ray sums at contact for EC emitters (Pb-203, Ba-133,
  Eu-152's EC branch).

### Conditional estimator: sum-peak-fed and x-ray-fed summing-in
`[low]`
- **Effect:** `CascadeMethod::Conditional` captures gamma-gamma sum-fed summing-in and
  vacancy x-ray summing-out, but not γ + coincident x-ray feeding into a window and it treats
  partners of 3+-gamma cascades as independent: ≈1–2% from FullRealization at contact
  (Co-57 136: 1.18 vs 1.20; Ba-133 356: 0.517 vs 0.494); ε²-small far away.
- **Implementation:** enumerate γ + x-ray-pair and triple-γ channels in
  `cascade_sum_pair_channels`, sample partners jointly, or route contact-geometry x-ray-heavy
  summers to FullRealization.
- **Shows up:** high-Z low-energy x-ray summers or sum-fed peaks at contact.

### IC-electron summing: residuals and an auto mode
`[low]`
- **Effect:** with `enable_ic_electrons` (default OFF, GEANT4-validated to <1% at 2 cm):
  (a) M+-shell and sub-15 keV conversion electrons are not deposited, under-counting
  summing-out on IC-heavy low-energy peaks at bare-source contact (Ba-133 P356 ≈ +14% vs
  GEANT4 there; ≈0 by 2 cm) — a crude full-energy M+ deposit over-corrected and was
  reverted; (b) nothing turns it on by proximity, and it costs ~5–8% whenever on.
- **Implementation:** (a) M-shell binding plus a partial-deposit model; (b) a tri-state
  `Off/On/Auto` (Auto within a few cm), re-tuning the cascade test bands, which are
  OFF-tuned.
- **Shows up:** (a) bare/thin sources at contact or near 4π; (b) close-source callers.

### L1/L2 conversion electrons use the L3 binding energy
`[low]`
- **Effect:** `build_radioactive_emissions()` and the IC-electron deposition take the L
  binding from `LFluorescenceData::l3_edge_keV` for L1, L2 and L3, so L1/L2 electrons are
  too energetic by up to ~4.8 keV (Np) or ~2.8 keV (U).
- **Implementation:** carry per-subshell L edges in `LFluorescenceData` (the EADL generator
  already reads them) and use the one matching each `ConversionShell`.
- **Shows up:** IC-electron summing or escape studies for strongly L-converted heavy-element
  transitions (Am-241, U/Pu chains).

### Coincident L x-ray sum peaks slightly low
`[low]`
- **Effect:** L x-ray singles match GEANT4 to ~0.3% (Am-241/Np, Pb-203/Tl), but coincident
  L sum peaks run low (Am-241 59.5+L ≈ 0.80, L+L ≈ 0.96 vs GEANT4), pointing at
  under-weighted upper-level feeding into the 59.5 keV level.
- **Implementation:** check the flow-conservation feeding in the level-path builder
  (`SandiaDecayCascade.cpp`, `feeding = out_flow − in_flow`) against Am-241's α branching to
  the Np 102.96 / 158 keV levels.
- **Shows up:** sum-peak magnitudes for high-Z L-x-ray emitters at close geometry (<0.25% on
  peak efficiency corrections).

### β⁺ annihilation in the cascade: residuals
`[low]`
- **Effect:** the cascade β⁺ path ranges the positron (allowed-spectrum KE, in-flight
  annihilation, escape from the source). Remaining: the source walk uses the electron
  (Møller) stopping term — `walk_in_source_geometry` takes no `is_positron`, although the
  Bhabha F⁺ term exists (`range_table_pos_`); no ortho-positronium 3γ; the Conditional
  estimator still annihilates at the vertex.
- **Implementation:** thread `is_positron` into `walk_in_source_geometry` (<1% effect, no
  cost); range the positron in Conditional.
- **Shows up:** β⁺ emitters (Na-22, F-18) at close geometry — the 511 / 1022 / γ+511 peaks.

### Internal-pair formation for verified high-energy E0 transitions
`[low — physics]`
- **Effect:** verified pure-E0 transitions get `p_gamma = 0`; above 2mₑc² the remainder may
  be internal pair formation, kept as an explicit `p_unmodeled` outcome, so N-16's
  6048.2 keV 0+ → 0+ transition no longer emits an impossible photon but its pair deposit
  is absent.
- **Implementation:** emit a back-to-back 511 pair plus the (E − 1022) keV kinetic energy,
  with the pair/conversion split from data. Do not promote GEANT4 multipole code 1 alone to
  an E0 assertion (Ac-232 373.3 keV and Rb-78 189.8 keV must stay photons).
- **Shows up:** verified E0 transitions above 1022 keV.

### High-coverage partial graphs keep small topology-free residuals
`[low — nuclear-data policy]`
- **Effect:** accepted graphs have at most 0.01 absolute residual occurrence (and below 99%
  coverage need complete RDM starts); residual transitions are kept as independent
  categorical transitions. Th-234 stays on its inferred DAG: its old RDM starts disagree
  with the newer 63.29/92.8 keV yields, and the record lacks GEANT4's 29.5 keV and 20.01 keV
  transitions, so 103.42 keV is a 0.0565/decay sink whose Pa L conversions are never
  generated.
- **Implementation:** SandiaDecay data (add both Th-234 transitions, re-derive the feeds, a
  reviewed levelFeed override); replace the independence approximation only with mapped
  topology or authoritative joint correlations.
- **Shows up:** unmatched transitions carrying <1% of a branch's occurrence sum.

---

## Source geometry

### Extended-source totals: a small coherent offset
`[low — fidelity]`
- **Effect:** with source-electron transport on and both sides at high precision, extended
  sources show sub-1% TOTAL structure that is statistically resolved: the uncored soil
  shell G-C reads +0.62% (z 2.1) at 238.6 keV, +0.74% (z 2.9) at 583.2, −0.63% (z −3.6) at
  2614.5; across the 16 cored + control rows χ²/n ≈ 3.2 with a positive mean (FEP is clean,
  χ²/n ≈ 1.1); config 12 totals sit +0.9 to +1.3% at 200–1173 keV. Cored and uncored cases
  agree, so it is the extended-source path, not cores.
- **Implementation:** unknown; candidates are the secondary-channel normalization above the
  pair threshold and the scattered-stream placement documented for config 8.
- **Shows up:** extended-source TOTAL at high statistics.

### Position bias ignores hollow-source voids
`[low — coverage]`
- **Effect:** `sample_biased_position` / `compute_auto_bias_params` use the outer
  `radius` / `half_dims` only, so position bias on a hollow source emits from the void with
  wrong weights. Unbiased sampling (the default for these shapes) is correct.
- **Implementation:** reject-and-reweight (V_solid/V_shell), annulus-aware sampling, or
  assert position bias off for hollow sources.
- **Shows up:** hollow sources large enough for auto position bias.

### Source-core validation gaps
`[low — validation coverage]`
- **Effect:** a partly filled core (leaving a real void) has no GEANT4 anchor; the G-E
  confinement control exists only for the sphere; nothing cored is anchored at ~50 cm (the
  far-field cone/two-stream path for a cored source).
- **Implementation:** a G-D-like half-filled case at near distance and G-E-style controls
  for the cylinder and box (a handful of 16M runs).
- **Shows up:** cores that leave a central void, or cored sources far away.

### A source overlapping the detector is simulated without complaint
`[low — API]`
- **Effect:** `EfficiencyCalculator` transports whatever source position and geometry it is given;
  a volumetric source cutting into the housing (or a point inside it) runs silently.
  `outermost_extent_radius()` only sizes the GDML world.
- **Implementation:** a `SourceGeometry::outer_solid(position)` checked with `convex_separation()`
  against `GeometryDescriptor::outer_envelope().solid()`. It needs an overlap depth, not the
  separation alone: a source at contact and one barely overlapping (by up to ~1e-4 of the size,
  near tangency) read the same. One way is to test the source solid shrunk by a tolerance. It must
  skip Marinelli sources, which wrap around the detector.
- **Shows up:** a host passing a mis-placed source.

### `add_core()` overflow is silently clamped
`[low — API]`
- **Effect:** core thicknesses summing past the cavity are clamped (`.cwiseMax(0.0)` in
  `rebuild_layer_dims`), silently dropping inner layers.
- **Implementation:** assert or throw when the accumulated thickness exceeds the cavity.

### Semi-infinite half-space (wall / floor / soil slab)
`[medium — coverage + estimator]`
- **Effect:** large-area surface or volumetric sources can be approximated today by a large,
  deep box, but uniform box sampling wastes nearly every history far off-axis without
  `lambda_lateral` position bias.
- **Implementation:** a half-space template with built-in lateral importance weighting,
  analytic lateral truncation, per-area/per-mass normalization, perhaps an angular-flux
  estimator. Air matters more here than anywhere.
- **Shows up:** in-situ / large-area measurements with standoff.

### More source shapes and activity distributions
`[low — coverage]`
- **Effect:** not expressible today: frustum / cone; bottle with neck and partial fill;
  several separated active regions in one run; graded activity (radial gradients, surface
  plate-out); activity in more than one layer (a trace source spread through a shield or
  attenuator — today one emitting layer sits between cores and shields). A thin disk /
  planar deposit already works as a thin cylinder (only a convenience constructor and a
  "massless deposit" flag are missing).
- **Implementation:** frustum needs a cone-intersection routine; bottle is a composite;
  multi-region and multi-layer activity need the source loop to sum over emitting regions;
  graded activity extends the sampler weighting (`PositionBiasConfig` /
  `DepthDistribution` hooks exist).
- **Shows up:** specialised containers; composite or graded sources; contaminated shielding.

---

## Detector geometry and GDML export

### GDML export emits no rotation for source volumes
`[medium — silent mismatch]`
- **Effect:** `write_gdml()` places source volumes with a position and an identity rotation,
  so a `CylinderSideOn` source or a rotated box is exported unrotated while the MC honours
  the rotation.
- **Implementation:** write the outermost source volume's rotation from `cyl_rotation_` /
  `rect_rotation_`; nested daughters share its frame.
- **Shows up:** validating a side-on cylinder or tilted box against GEANT4.

### Dead layer is not exported to GDML
`[low — validation coverage]`
- **Effect:** `write_gdml()` writes the whole crystal as `active_crystal` and only notes the
  dead layer in a comment, so a dead-layer geometry is compared against a different solid.
  Configs 25/26 carry no dead layer, and the bulletized active fillet has no GEANT4 anchor.
- **Implementation:** export the crystal as a non-scoring volume with the active region as a
  daughter (same material, so not a subtraction); the harness already scores
  `active_crystal`.
- **Shows up:** HPGe with realistic 0.5–1 mm dead layers, especially below ~100 keV.

### Side & back incidence: no near-field model behind the face plane
`[medium — accuracy close beside/behind]`
- **Effect:** `back_coverage` scans the far field only; close beside or behind the detector the
  kernel alone carries the distance dependence, with the `NearFieldUnmodeled` envelope; and
  behind a thick back, the removal-only kernel misses full-energy in-scatter (shadow-gated, so
  flagged rather than wrong, but unanswered).
- **Implementation:** a back near-field table on (E, cos theta < 0, phi, clearance s/a), clearance
  measured past `outer_envelope()` so no node falls inside the detector, with an ln N = 0 anchor
  row. In the shadow regime (s < 0.3) even back nodes do not help - the in-scatter ratio varies too
  steeply between them (measured 20-70% off) - so that needs either much denser nodes there or a
  scatter term in the kernel. `closed_loop` probe families and transfer-mode back anchors could
  follow.
- **Shows up:** sources within a few detector widths of the side wall or the back.

---

## Electron / CSDA physics

### Retire the empirical skin-escape gate with a Mott / Goudsmit–Saunderson MSC
`[low — robustness; research]`
- **Effect:** the regime-aware exit-state gate is the validated default (mean |err| at
  662 keV: gate 0.56% < `tails` 0.62% < `gs` 0.92%; DESIGN.md has the per-config table). The
  first-principles variants over-diffuse multi-wall boxes (`tails`) or under-scatter ~30%
  (`gs`, screened Rutherford without Mott and nuclear form-factor corrections).
- **Implementation:** replace the screened-Rutherford cross section with Mott, or sample a
  tabulated Goudsmit–Saunderson angular distribution by inversion (no Gaussian+tail
  double count); then re-validate ~1% across Z without the gate. Run time is ~free; the
  cost is developer time and Mott data.
- **Shows up:** shielded/extended box sources (configs 8, 11, 12, 20–24).

### Electron tables for Np–Cf are uranium's
`[low — coverage]`
- **Effect:** ESTAR, ICRU-49 I and EPQ Seltzer–Berger stop at Z = 92, so Z 93–98 reuse
  uranium's (DESIGN.md table: collision −2.1 .. +2.3%, radiative −1.9% (Np) .. −8.9% (Cf),
  CSDA range −1.7 .. +2.9% above 50 keV). The efficiency effect is unmeasured.
- **Implementation:** add ESTAR 93–98 as a production source and a bremsstrahlung source
  above Z = 92 (EEDL2023 MF=26 MT=527 would be a new lineage). Cheaper for the radiative
  part: rescale uranium's by Z(Z+1)/A (error −2.5 .. +3.9%).
- **Shows up:** thick Pu/Am/Cf sources or shields at MeV energies.

### In-crystal positron annihilation: residual approximations
`[low]`
- **Effect:** the PP positron is CSDA-tracked with Berger–Seltzer F⁺ stopping, annihilates
  at its endpoint or in flight, and the pair gets the Koch–Motz opening angle (single/double
  escape peaks now on GEANT4: −1.3% / −0.5%). Remaining: the in-flight pair uses an
  isotropic energy split rather than exact two-body kinematics (a candidate for the small
  single-escape offset); no ortho-positronium 3γ; the in-flight cross section is duplicated
  (`positron_annih_xsec_per_electron` in `ElectronCsda.cpp`, `annih_xsec_per_electron` in
  `EfficiencyCalculator.cpp`).
- **Implementation:** exact boost kinematics; one physics home for the cross section.
- **Shows up:** escape-peak / annihilation-continuum shape above 1.5 MeV, thin crystals.

---

## Cross-section data

### Atomic weights through uranium are xraylib's two-decimal table
`[low — fidelity]`
- **Effect:** `g_atomic_weights` holds xraylib 4.2.1's rounded, partly outdated values: H
  1.01 (+0.20% vs IUPAC), Au +0.12%, Ti +0.07%, Mg +0.06%, Ge −0.06%, Na +0.04%, Al
  −0.04%, and Tm 167.27, which is erbium's (−0.98%). Atom densities move by the same
  fraction the other way (hydrogen in water −0.2%, electron density −0.04%). The GDML export
  uses IUPAC masses, so CeeLo and GEANT4 already differ by these amounts.
- **Implementation:** take IUPAC 2021 abridged weights for Z ≤ 92 as a locked source in
  `generate_element_support.py` (keep the actinide conventions below), re-run the migration
  validators and spot-check the gates.
- **Shows up:** hydrogenous and HPGe materials at 0.05–0.2%; thulium at 1%.

### Actinide atomic weights are one conventional mass per element
`[low — coverage]`
- **Effect:** above uranium CeeLo uses xraylib's masses (Np 237, Pu 239.1, Am 243, Cm 247,
  Bk 249, Cf 251). Pu and Bk match the isotope in hand; a material that is really Am-241
  gets a μ 0.8% low (Cm-244 1.2%, Cf-252 0.4%).
- **Implementation:** let `Material` take a per-component atomic mass or nuclide
  composition (InterSpec knows the isotopics), keeping xraylib's as the default.
- **Shows up:** per-gram quantities for Am/Cm/Cf sources or shields.

### Rayleigh anomalous dispersion (f′, f″)
`[low]`
- **Effect:** the coherent angular model uses non-dispersive F(x,Z), overestimating by
  several to tens of percent in narrow windows above high-Z edges. EPICS2023's *integrated*
  coherent cross section is 1–34% below XCOM's for Z 92–98, mostly just below K edges — the
  anomalous-scattering pattern — so only the angular sampling may lack it.
- **Implementation:** f′(E) + if″(E) (xraylib Fi/Fii or RTAB) in the Rayleigh angular model,
  after checking the integrated term.
- **Shows up:** sub-30 keV high-Z work; Rayleigh-peak shape.

---

## Error estimation

### No systematic-error budget
`[low]`
- **Effect:** reported σ is statistical only, while documented model biases (LaBr3 ≥2 MeV,
  cfg 8 at 59 keV, skin escape) and data differences (cfg 12 at 59 keV, cfg 27/28 at
  60 keV) exceed σ_stat in their regimes, so z = Δε/√(σ_MC²+σ_G4²) is spuriously large
  there. Cfg 12 at 59 keV (FEP −4.0%, total −3.2% against tolerances 2.1% / 1.9%) passes
  the gate only on its 3σ allowance — a known offset carried as noise.
- **Implementation:** document (optionally expose) an energy/config-dependent σ_sys; at
  minimum widen or `SKIP` cfg 12 at 59 keV with the reason attached, as cfg 7 at 100 keV is.
- **Shows up:** MC-vs-GEANT4 z-scores in the known-bias regimes.

---

## Minor

- Rename the `MCDET_NO_{STRAGGLE,B2,TAILS}` env knobs (`ElectronCsda.cpp`) to `CEELO_`.
- `PhotonTransport.h` still mentions `diag_interactions`, now a dev-only study tool.
- `ElectronCsda.h`'s `kSB_min_energy_keV` comment says the brems tables start at 200 keV;
  they start at 10 keV.
- Doppler broadening omits the F(p_z) profile-reweighting factor (few-% tail) and the
  relaxation cascade from the Compton vacancy — Compton-edge **shape** only; never cite it
  for FEP/total.
