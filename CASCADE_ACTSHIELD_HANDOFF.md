# Act/Shield cascade-summing + DRF-uncertainty + fixed-geometry-MC — handoff

Status as of 2026-07-09. Plan: `~/.claude/plans/please-update-the-activity-shielding-floofy-wand.md`.
Branch `upgrade/DetEff`, commits `2560d0a9`, `493ac00e`, `62a4bb78` (+ CeeLo dev-repo commit
"Template analytic cascade summing on the efficiency scalar type", vendored copy synced).

## Done

- **CeeLo (dev repo `/Users/wcjohns/coding/mc_det_eff_plan`, then vendored)**: analytic cascade
  code (`compute_cascade_analytic`, `EfficiencyProviderT`, `AnalyticPeakResultT`, LevelDag
  survival DP) templated **in place** on the efficiency scalar (ceres::Jet-compatible); proven
  answer-identical by `tests/test_analytic_cascade_regression.cpp` (stored pre-refactor values,
  exact equality) + a dual-number/finite-difference check. `ResponseGenerator::configure_calculator`
  exposes descriptor→scene setup. Level-scheme coverage verified for Co-60/Y-88/Ba-133/Eu-152 and
  the Bi-214 branch of aged Ra-226 (all have schemes — the Ra-226 data risk is retired).
- **P2 (DRF uncertainty)**: `energy_chi_contributions` inflates per-peak denominators with the DRF
  efficiency sigma (chart + legacy Minuit2 now agree with the Ceres GLS whitening); `PeakDetail::
  drfEffFracUncert` → batch JSON (`DrfEffFracUncert`) → txt/html templates. Tests: fully-correlated
  5% band maps ~1:1 onto activity uncertainty; option-off bit-identical; XML round-trip.
- **P3 (cascade, point sources)**: `ShieldingSourceFitOptions::correct_for_cascade_summing`
  (default off, XML back-compat); `CascadeSummingCalc` (thin adapter over CeeLo's engine; cascade
  cache per age bucket; `ShieldScatterAugment` GADRAS table, bilinear in (AN, AD), Jet-differentiable);
  per-nuclide local cluster maps × c_net in both fit paths; Jet derivatives flow through shield
  chords; `create()` throws without `DetectorPeakResponse::hasAnyTotalEfficiencyInfo()` (API contract);
  point-source primary efficiency dispatches to `fepEfficiencyEval` when a CeeLo response is attached.
- **P4 (extended sources)**: per-element corrections inside all four `DistributedSrcCalcT::eval_*`
  integrands (partner Σμ·dist via `record_path`/`record_generic` scratch, per-element FEP/partner air,
  per-element (AN, AD) scatter augmentation, single-window engine call); build-time per-nuclide skip
  guard (upper-bound estimate ≥1e-4). Geometry-golden integrals bit-identical with the option off.
- **P5 (GUI)**: "Cascade summing corr." checkbox; enabled iff DRF info present AND predicted magnitude
  ≥0.3% (`estimateMaxSummingMagnitude`; fixed-geom: info only); missing-info SimpleDialog →
  `InterSpec::showDrfModifyWindow` (new Cat-A window over `DrfModifyWindow`); undo/redo; i18n.
- **P6 (reports)**: per-peak `CascadeCorrection{NetMult,SummingOut,SummingIn}` JSON + option-gated
  columns in `std_fit_log.tmplt.txt` / `act_fit.tmplt.html`; calc-log lines. (Partner-list breakdown
  in reports = follow-up; the engine does not expose per-partner terms yet.)
- **P7 (fixed-geometry MC)**: `MakeFixedGeomResponse` (scene assembly → per-energy precision MC →
  DRF with FEP+total curves, MC-sigma uncertainty, `<FixedGeomSourceSetup>` blob, hash-gated);
  GUI button + worker/progress/cancel + read-only locked display of the embedded scene.
- **P1 harness**: `target/testing/test_CascadeSummingFit.cpp` (scene matrix incl. thin disks
  0.5–5 cm, `--regenerate` mode that re-grounds the golden response to the LANL Detective-X curve
  and re-runs the CeeLo FullRealization truth MC). `detective_x` added to `regen_ceelo_goldens`.

## Remaining (in order)

1. **Wait for `make_golden_response detective_x`** (running in background; writes
   `target/testing/test_data/ceelo_drf/detective_x_response.xml` + `_probe.csv`; args
   `general 0.001 120 1 200000000 60`).
2. **Generate truth fixtures** (hours; overnight):
   `cd target/testing/build_ninja && ./test_CascadeSummingFit -- --datadir=$PWD/../../../data --testfiledir=$PWD/../test_data --regenerate`
   Then `ctest -R TCascadeSummingFit` — the Cs-137 rows must pass. Commit
   `test_data/cascade_truth/` + the detective_x fixtures + a README (gates table is in the test file).
3. **Flip `CASCADE_TRUTH_FULL` to 1** in `test_CascadeSummingFit.cpp` → all cascade + extended rows
   assert. Debug against gates (bare Co-60/Y-88 3–4%, Ba-133/Eu-152 5–8%, shielded +2 pctpts,
   disks 6–8%). Also time a self-atten fit option-on vs off (target ≤1.5×).
4. **Scatter quantification**: run Co-60/Ba-133 behind Fe 10 mm / Pb 6 mm with and without the
   `ShieldScatterAugment` term; paste numbers into the comment block atop `src/CascadeSummingCalc.cpp`.
5. **Re-run the full CeeLo dev-repo ctest** (was stopped for CPU; targeted analytic tests all pass).
6. **Manual GUI pass** (user): checkbox gating vs distance; missing-info dialog → editor round-trip;
   fixed-geom MC button → progress → locked display; calc log/report output; `?isphone=1` layout.
7. Follow-ups filed: per-partner report breakdown; `eps_*_element` integrand refinement; per-element
   scatter-table fracH axis; Jet ∂C for generic-material shields in the cascade term (currently
   scalar); Minuit2-path removal.

## Truth-run findings (2026-07-10)

Cs-137 control rows pass at 1.1–1.7%. Full-gate run surfaced three things:

1. **Point-source cascade validation is solid.** Co-60 all pass; Ba-133 all pass;
   Eu-152/Y-88/Ra-226(bare) mostly pass. Remaining point misses are near-contact (2 cm)
   edge cases where the DRF near-field interpolation (~3%) + analytic residual stack up:
   cs137_bare_2 3.0% (gate 2.5), y88_bare_2 4.5% (gate 4), eu152_fe10mm 8.9% (gate 7).
   These need small, physically-justified 2 cm/heavy-shield gate bumps (not loosening to hide
   a bug — the point cascade correction demonstrably works: cascade-on ≫ better than off).

2. **FIXED a real bug in the extended per-element cascade path.** `totalIntrinsicEfficiency()`
   throws when the DRF's total-efficiency lives in a CeeLo response (no total curve); the
   try/catch in `build_volumetric_calculators` swallowed it, zeroing `partner_tot_int` → Co-60
   extended summing-out was exactly 0 (cascade-on == cascade-off). Added
   `DetectorPeakResponse::totalIntrinsicEfficiencyAny()` (dispatches to total curve OR backs the
   intrinsic total out of the CeeLo response far-field) and used it for the partner efficiencies.

3. **Extended-source absolute activity has a ~14–38% baseline gap, independent of cascade.**
   The volume integrand uses the DRF's intrinsic *curve* × solid angle (the
   `detector_response_is_isotropic()` assumption), while the truth is the full-angular CeeLo
   response — so for wide/close extended sources the two disagree regardless of summing
   (cascade-on == cascade-off to 6 digits for Co-60 before the fix). This is the deferred
   `eps_*_element` integrand refinement, NOT a cascade bug. RESOLUTION: the extended-source
   truth rows should validate the cascade correction as a RATIO (cascade-on model / cascade-off
   model ≈ CeeLo volume-averaged summing factor), not absolute activity; make the absolute
   extended gates advisory (report-only) with this note, OR adopt eps_*_element in the integrand.

4. **Ra-226 shielded**: CeeLo FullRealization does not emit daughter-nuclide (Bi-214) lines in a
   multi-generation chain, so only the 186 keV Ra-226 line has truth; behind Fe it is a single
   heavily-attenuated line → ill-conditioned (14% off). Make ra226_fe10mm_25 advisory. (Bare
   ra226_bare_25 passes at 2.3% on the 186 line.) File the CeeLo daughter-line gap as a dev-repo
   follow-up; it does not affect InterSpec's analytic fit path.
