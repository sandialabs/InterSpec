# CeeLo — Open Items / TODO

The single register of known shortcomings, deferred features and fidelity/speed work.
**When an item is resolved, delete it in the same commit that fixes it** (and update
DESIGN.md if it described the gap).

Each item gives its **Effect** (what is wrong or what it buys), **Implementation** (what it
would take) and **Shows up** (the regime where it matters), plus a severity tag relative to
the ~1%-vs-GEANT4 goal. Areas are listed in priority order — the full-spectrum
approximation first on this branch — and items within an area are sorted the same way;
"Next up" is the short cross-area ranking. Items measured to be sub-statistical are kept
for awareness. Full-spectrum payloads are regenerated and accepted with
`tools/full_spectrum/README.md`; "the last bake" below means the most recent bake of the
current representation (where a gate has no documented bar, its bound is that bake's value
with headroom).

## Next up

Full-spectrum:
1. [Move the detector validation bank off the training grid](#the-detector-validation-bank-sits-on-training-points)
   — its full-energy bar is below its own noise, so the gate cannot be trusted yet.
2. [A local-linear yield evaluator](#the-full-energy-atom-between-nodes) — the full-energy
   gate fails on every detector; a format change, so re-bake all detectors after it.
3. [Solve the continuum-shape nodes](#the-detector-bakes-miss-arrival-gates) — the arrival
   gates.
4. [The campaign-vs-`EfficiencyCalculator` low-energy FEP offset](#the-detector-campaign-and-efficiencycalculator-disagree-below-150-kev)
   — cheap, and the composition's absolute scale depends on it.
5. [Commit the Nov 3 source-escape and charged-particle payloads](#commit-the-nov-3-source-escape-and-charged-particle-payloads),
   then [a payload loader](#no-payload-loader-full_spectrum_nuclides-ignores-the-source-payload).
6. [The source → environment → detector composition layer](#no-source--environment--detector-composition-layer).
7. [Near-field material and detector characterization](#near-field-material-where-a-fold-stops-being-valid)
   and [environment accuracy targets](#environment-accuracy-targets-and-the-adaptive-round-cost)
   — the next two program steps.
8. [Real Detective X / Micro-Detective payloads](#no-fitted-payload-for-a-real-field-hpge), then
   [NaI](#no-nai-detector-payload).

Shared engine (as on main): [the GEANT4 gate in ctest](#compare_validationpy-is-not-wired-into-ctest),
[load-robust precision tests](#precision-targeted-mc-tests-are-wall-clock-capped),
[the 0.75 keV references](#geant4-references-are-still-scored-at-the-old-15-kev-fep-window),
[the Na-22 511 keV drift](#na-22-511-kev-drifted-from-geant4-after-the-positron-model-change),
[decay-data level schemes](#rejected-branch-pair-links-cannot-all-be-honored),
[GDML rotation](#gdml-export-emits-no-rotation-for-source-volumes),
[air](#air-between-source-and-detector).

---

## Full-spectrum: detector response

### The detector validation bank sits on training points
`[medium — validation]`
- **Effect:** 56 of the 88 validation cells share their (E, cos θ) with `dense_energy`
  training cells (only the seed differs), so bandwidth selection and the validation gates
  mostly score on-node reproduction, not interpolation. The full-energy p95 bar (1.5%) is
  also below the bank's own noise: above 3 MeV a cell holds 2–60 thousand full-energy counts
  (σ 0.4–2.3%), so a perfect model would read `hpge_coax`'s p95 as 1.71% (90% range
  1.22–2.40%) and fail the bar 73% of the time.
- **Implementation:** move the validation energies and cosines off every training grid (the
  2% rule of `detector_bank_energies.py`, as the audit and arrival banks are), give cells
  above 3 MeV enough events for σ ≤ 0.4% (~6M at 11 MeV; ~1 CPU-hour per detector), and
  re-select the bandwidths on the new bank.
- **Shows up:** reading the validation gates as between-node accuracy; every bake's
  full-energy p95 gate.

### The full-energy atom between nodes
`[medium — accuracy]`
- **Effect:** with the yield nodes solved to reproduce their training cells, `hpge_coax`'s
  full-energy atom reads 0.28% / 2.29% (median / p95; bars 0.5 / 1.5%) on validation and
  0.46 / 2.89% and 0.35 / 2.10% on the audits; `detective_x_fit` 2.28%, `detective_x` 1.65%
  p95 — every detector fails. The misses are between nodes, where a kernel average of node
  coefficients cannot follow a steep gradient: every 11 MeV cell reads low (−1.5 to −4.9%,
  mean −2.1%, z ≈ −4.7); near the axis (|cos θ| 0.93–0.97) −3 to −4.5% at 0.5–10 MeV (the
  yield changes steeply between θ = 0, a coordinate boundary, and the next node); and on the
  axis between energies, `detective_x_fit` reads −1.32% (z −23) at 344.28 keV and +0.59% at
  356.01 keV, which holds the review page's Eu-152 normalizing peak at 0.98 of Monte Carlo.
- **Implementation:** a local-linear yield evaluator (the source-escape secondaries already
  use one in C++): a kernel-weighted affine fit at the query follows a gradient across a gap
  and has no boundary bias at θ = 0. Fit the gradient in (log E, θ) only, with the kernel
  over all three coordinates (the passive coordinate is a function of energy); the node solve
  then applies to the local-linear smoother. A format change: regenerate every payload (the
  banks stay).
- **Shows up:** full-energy efficiencies near axial incidence off the training energies, or
  above 3 MeV; every bake's full-energy gate. The first `hpge_coax` bake with spectrum
  fidelity (Nov 3: Detailed walk, Compton K vacancies, impact ionization) reads p95 2.93%.

### The detector bakes miss arrival gates
`[medium — accuracy]`
- **Effect:** the arrival bank (15 held-out energies × iso / front / back / side) is the
  response a scattered field is folded with. `hpge_coax` misses one gate (front 40–100 keV
  detected-yield median 0.163% against 0.15%; the rate surface's). `detective_x_fit` misses
  five: detected-yield medians at 40–100 keV front and side, and valley-density medians
  3.3–3.4% at 400–1000 keV, where its 340–640 keV cells read +2.4 to +6.2% (z 6–19) in every
  class while the full-energy atom beside them is within 0.9% — the full-energy correction
  moved yield into the continua, whose shapes are still plain kernel averages (bandwidth
  ~9% in energy). The valley misses are deferred until real spectra say they matter. The
  first `hpge_coax` bake with spectrum fidelity (Nov 3) misses four, all at 40–100 keV:
  detected-yield median back 0.169% and side 0.177% (front passes), back valley-density max
  11.2% against 10%, side detected-yield excess max 0.86% against 0.5%; its bounds were set
  on fidelity-off truth.
- **Implementation:** solve the continuum shape nodes as the yield nodes are, keeping the
  coefficients nonnegative (log/softmax parameterization or projected steps), and re-score.
  At the same time audit the `multiple_compton` block boundary at the free-electron Compton
  edge (a genuine step for `single_compton` but not necessarily for multiple scatter; ≤3.2σ
  at 5 MeV, never audited after resolution folding): if a jump shows, constrain the two
  blocks to match in value (and slope) there.
- **Shows up:** folding a scattered field below 100 keV or into the 0.25–1 E_line valley;
  every bake's arrival gate.

### The detector campaign and `EfficiencyCalculator` disagree below 150 keV
`[medium — consistency]`
- **Effect:** for `hpge_coax` head-on the payload reproduces its campaign truth (59.54 keV:
  −0.12% full energy), but `predict` at 10 m reads the full-energy peak +1.02% (z 26) at
  59.5 keV and +0.45% (z 19) at 122 keV above a 4M-event `EfficiencyCalculator` run of the
  same geometry, while detected totals agree to 0.1% and 662–2614 keV agree within 0.2%;
  `detective_x_fit` reads +1.4 to +1.7% at 59.5 keV. Both sides transport fluorescence and
  bremsstrahlung; the campaign counts a deposit within 1e-6 keV of the incident energy as
  full energy, `EfficiencyCalculator` anything within 0.75 keV — which would push the other
  way. Yet the review page's 59.5 keV peak fit to FS and MC spectra of Am-241 on that crystal
  gives FS/MC = 1.001.
- **Implementation:** tally both definitions in one run of each code at 59.5 and 122 keV and
  find the histories one counts and the other does not; then why the page's peak fit does
  not see it.
- **Shows up:** functional-model vs Monte Carlo full-energy efficiency below ~150 keV.

### No fitted payload for a real field HPGe
`[medium — coverage]`
- **Effect:** the one gate-checked HPGe, `hpge_coax`, is a generic 8 × 8 cm coax, and a
  payload covers exactly one detector. `micro_detective.geometry` (the InterSpec DRF "ORTEC
  Micro-Detective (CeeLo geometry, LANL 25 cm curve)", landmark 62.2 keV) has never been
  baked; `detective_x.geometry` is a **placeholder** (a Detective-X-like crystal in the
  Micro-Detective housing). The latest Detective X bakes (`detective_x`, `detective_x_fit`)
  fail only the full-energy p95 gate and arrival gates (above), now that the low-energy banks
  follow each detector's landmark; the bandwidths were still selected on `hpge_coax`.
- **Implementation:** replace `detective_x.geometry` with the real geometry (`drf_geometry.py`
  from a DRF that carries one; commit the geometry files the bakes use), then bake both
  (`DETECTOR=micro_detective`, `DETECTOR=detective_x`; ~1–1.5 h each on a 10-core Mac
  including the 25–60 min finite-distance step), re-checking bandwidths per detector. After
  the three items above.
- **Shows up:** any Micro-Detective or Detective X spectrum or FEP efficiency from the model.

### Point-source detected totals on `detective_x_fit` read 0.1–0.9% low
`[low — accuracy]`
- **Effect:** `predict` head-on reads the detected total −0.48 / −0.64 / −0.83 / −0.51% at
  10 cm / 30 cm / 1 m / 3 m (661.657 keV), −0.15 to −0.38% at 2614.511 keV and ~−0.85% at
  344/356 keV at 1 m against 4M-event `EfficiencyCalculator` runs (z −3 to −29), while its
  full-energy peak agrees within 0.3%; `hpge_coax` reads +0.03 to +0.61%.
- **Implementation:** run the ladder to 10 m on `detective_x_fit` to separate the payload's
  detected yield (rate surface) from the detected depth, then fix whichever it is.
- **Shows up:** total-count or continuum fits with a Detective X payload at 10 cm–3 m.

### No NaI detector payload
`[medium — coverage]`
- **Effect:** every fitted detector payload is HPGe; `nai3x3_payload()` (the default) and
  `coaxial_hpge_payload()` are hand-set pilots. The environment's detector-fold reports need
  a NaI response too. `DETECTOR=nai3x3` runs, but a smoke bake stops in the fit ("no
  populated nodes for annihilation_511" — perhaps only smoke statistics), and the recipe is
  HPGe-tuned (the `lowE` bank, the K-line components and their support limits were chosen on
  Ge K lines; NaI has I K x-rays at 28–33 keV and a far larger escape fraction).
- **Implementation:** run the NaI bake at production statistics, fix what the fit needs,
  re-check the bank energies against NaI's landmark (`ceelo_rate_surface_anchor --info`).
- **Shows up:** any NaI spectrum from the functional model.

### `predict` uses the far-field continuum shape for a nearby source
`[medium — accuracy]`
- **Effect:** `FiniteDistanceDepth` corrects the yields of the tracking atoms and the detected
  total for a point source's distance, but the continua keep the parallel-beam shape the payload
  was fitted to. Rays from a near source enter the face obliquely and leave through the side
  after a short chord, so more of the continuum sits at low energy: on `detective_x_fit` at 10 cm,
  CeeLo's own MC puts 8–13% more of its detected counts below 100 keV than the parallel beam at
  cos θ = 1 (662 and 1332 keV, z +3 to +8) and 2.5% less in the full-energy peak; at 122 keV,
  +34% at 0–10 keV and −9% at 60–120 keV. `predict` inherits exactly that: in the
  CeeLo-vs-GEANT4 study (`scratch/fullspec_g4`, 10 cm) FS / MC reads 0.88–0.93 below 100 keV
  and +1.3% from the backscatter shelf up (662 keV bare), and at 122 keV FS / MC equals
  beam / point within 1% in every region, while the payload reproduces its own training cell at
  661.657 keV head-on within 0.4% everywhere.
- **Implementation:** have the bake's `finite_distance` step tally the deposited spectrum too
  (it already runs point-source MC at 21 energies × 13 distances) and fit a smooth
  distance-dependent shape correction in the continua's landmark coordinates (a low-order
  multiplicative tilt per block, its amplitude following the virtual disk's distance law),
  applied in `predict` with the depth factors.
- **Shows up:** any source within ~1 m, bare or shielded: the continuum below ~150 keV and
  around the Compton edge.

### Finite-distance depths are measured on the axis only
`[low — coverage]`
- **Effect:** `predict` applies the on-axis virtual-disk depths for cos θ ≥ 0.5 and the
  far-field acceptance below. `hpge_coax` at 662 keV, 60°: full energy −0.72 / −0.45 /
  +0.03%, detected −0.21 / −0.29 / −0.11% at 30 cm / 1 m / 3 m; side-on the far field reads
  −1.7% at 30 cm; backward incidence is unmeasured.
- **Implementation:** measure depths at a few incidence angles (cos θ = 1, 0.5, 0, −0.5, −1),
  each with its own distance convention (a side-on source is nearer the crystal than d, so
  the depth may be negative), and interpolate in cos θ.
- **Shows up:** `predict` for a nearby source off the detector axis.

### Box crystals are rejected by the detector anchor
`[low — coverage]`
- **Effect:** `set_rate_surface_anchor_geometry` throws for a box crystal, so no payload can be
  baked for CZT, CdTe or GADRAS's box-shaped HPGe descriptions (the MC handles boxes). For
  GADRAS's Detective X (5.87 × 5.87 × 4.59 cm) an area-equivalent cylinder (R 3.312 cm)
  matches the box's MC peak efficiency within 0.2% at 59.5–122 keV and reads 2–3% high above
  300 keV.
- **Implementation:** box face / side / back projected areas and entry stacks in the anchor,
  and a box envelope in `predict`'s acceptance.
- **Shows up:** box-shaped detectors.

### Detector fits sit at the truth-bank noise floor
`[low — accuracy; needs MC]`
- **Effect:** the detected rate is at 1.0–2.4× the truth cells' own binomial σ (52.5–120 keV
  1.0×, 130–2000 keV 1.3×, >2 MeV 1.9×, 42.5–45 keV 2.4×); the charged-particle payload is
  similar. Refitting banks of this size mostly fits noise: the binding constraint is
  training statistics, not capacity, regularization or coordinates (DESIGN.md, Previously
  tried).
- **Implementation:** more statistics on the 42.5–45 keV cells (σ 0.13–0.31% against
  0.05–0.15% elsewhere) — the only group with real headroom.
- **Shows up:** pushing the detected rate below ~0.2%, or ~0.5% at 42.5–45 keV.
- **Continuum shapes too:** near a training energy the shape field returns that cell's own
  coefficients, MC noise included. At 661.657 keV head-on the payload's residual against a
  32M-event spectrum correlates 0.8 with the training cell's (500k events) own noise, 1.6% rms
  at 10 keV scale, which shows as ripples in a raw spectrum. Remove them with a roughness
  (second-difference) penalty in the per-cell shape fit, chosen against held-out high-statistics
  cells, and more events per training cell (16× takes the ripple to ~0.4%; bake CPU is not the
  constraint); gate it with the detrended residual against a high-statistics reference. Where a
  continuum is continuous across a block boundary (multiple Compton at the free-electron edge,
  the pair blocks), tie the blocks in value and slope so the boundary cannot step.

### GEANT4 detector anchors lack the K-line labels
`[low — validation]`
- **Effect:** the harness's component table (`RunAction.cc`, `component_names[]`) has the
  original twelve labels, so a GEANT4 anchor lumps the K-escape and absorbed-K lines back into
  `photoelectric_escape` (and, after a Compton scatter, into the Compton continua, where the
  CeeLo campaign now counts them as atoms too); the matched-physics anchors (800 / 2614.5 / 5000 keV) also predate
  the crystal-walk straggling changes, and their +1.44% (z 4.3) line residual at 5 MeV has
  not been re-measured.
- **Implementation:** add the labels and their classification (the crystal fluorescence-escape
  flag and the dead-layer-absorption case) to the harness, rebuild, regenerate the anchors at
  42.5–120 keV and the three high-energy points — together with the physics-list A/B
  (Validation area).
- **Shows up:** component-resolved GEANT4 anchor checks below ~120 keV.

---

## Full-spectrum: environment

### No source → environment → detector composition layer
`[high — feature gap]`
- **Effect:** `FastSpectrumModel::predict` composes source escape → one detector incidence →
  resolution **in vacuum**; the environment payload (`src/io/EnvFitPayload.h`) returns
  incident φ(E, μ) at a detector sphere and folds nothing. A spectrum with air or ground
  scatter cannot be evaluated end to end. The environment's Rayleigh/direct ownership rule (a
  Rayleigh-on payload pairs with the `AnalyticCollimated` direct term) is a convention, not a
  type (the query takes an unlabelled `direct_fluence_per_cm2`), and the interface keeps one
  arrival cosine and assumes isotropic point-like leakage without stating or validating it.
- **Implementation:** one scene-level API: the direct term with analytic air attenuation; per
  escaping line/continuum bin, the environment payload's scattered φ(E, μ); both folded
  through the detector response with `ArrivalResponseTable` (it exists and is Jet-safe; its
  accuracy is in DESIGN.md), keeping atoms as atoms until resolution. Enforce the
  Rayleigh/direct rule in types, reject non-isotropic leakage and out-of-domain queries,
  return direct / air / ground / detector parts separately, and test the absolute scale
  against `EfficiencyCalculator`. It can be built now against `hpge_coax`.
- **Shows up:** any spectrum that includes air or ground scatter.

### Near-field material: where a fold stops being valid
`[medium — accuracy, coverage]`
- **Effect:** folding the analytic single-scatter flux of a plate 0.5–10 cm behind the crystal
  through the exact parallel-beam response reads 0.48–0.83 of a full MC with the plate in the
  detector geometry (z −3 to −20; shape TV 38–100%; 122–1332 keV, 30–300 cm): crystal → plate
  → crystal summing moves counts that no incident-flux fold represents, so such material
  belongs in the detector geometry. Where a fold becomes valid is unmeasured (plates 20–40 cm
  behind add 0.1–2% of the FEP count, below what two independent 4M-history runs resolve),
  and how far near-field clutter is degenerate with shielding and activity in a fit is not
  yet quantified.
- **Implementation:** correlated sampling or plate-tagged events to find the fold-validity
  distance; carry detector-attached material (housings, benches, collimators) in the
  payload's geometry; quantify the degeneracy with synthetic fits.
- **Shows up:** spectra with structures within tens of centimetres of the detector.

### Environment accuracy targets and the adaptive-round cost
`[medium — bake cost; decision]`
- **Effect:** after the field/soil baseline (384 nodes × 4 blocks, 13.3 GB) no node meets its
  noise target (worst uncertainty ratio median ~6.5, i.e. ~40× the histories at the median
  node; 115 of 256 training nodes need the 244-block cap or more). The capped round-2
  projection is ~67,500 more blocks (17.7 G histories, ~80 h of transport, ~13 GB more
  archive, at most ~55 GB more; ~27 rounds at 64 nodes per pass). 5 of 320 nodes hold 52% of
  the loss, and component maxima of 60–200% yield error are far above any node's MC σ, so a
  model part remains to be seen. The targets were never derived from what a field fit needs.
- **Implementation:** set the targets from synthetic field fits first (planned: smooth TV
  ~0.3%, MC noise 0.1%, environment model error carried into fits as a covariance), then run
  the rounds.
- **Shows up:** planning the environment campaign.

### Environment Doppler broadening is off down to 20 keV
`[low — shape, low energy; decide before the rounds]`
- **Effect:** the production campaign covers E0 = 20–3000 keV with `doppler = false` pinned
  by `publication.py`. That is harmless at hundreds of keV but not near 20–40 keV, where the
  once-scattered continuum sits within ~1–3 keV of E0 and Doppler smearing (~0.3–1 keV) is
  comparable. At 59.5 and 122 keV (1 m / 1 m, 2M histories) folded through an HPGe response,
  on vs off changes scattered counts ≤3% per 4 keV bin, within noise; 20–40 keV is
  unmeasured.
- **Implementation:** expose `EnvScatterConfig::doppler` in `env_kernel_bake`, A/B a 30 keV
  node through HPGe resolution, and turn it on if it matters (a new campaign identity, so
  before the rounds).
- **Shows up:** environment-scatter shapes for sources below ~60 keV.

### Environment-scatter publication is blocked
`[medium — publication blocker]`
- **Effect:** `promote-candidate` needs `CEEFIT-DETECTOR-RESPONSE-v1` NPZ matrices (HPGe and
  NaI) that nothing produces (`mu-fold-report`, `detector-fold-report`), a
  `CEEFIT-GEANT4-ANCHOR-MANIFEST-v1` that nothing builds (`geant4-anchor-report`), and the
  full-statistics registered-E0 gate (`env_registered_e0_gate`, 25 min–2 h), never run. The
  bake's own campaign would also be refused: `bake/env_scatter.sh` designs on 20–3000 keV,
  0.25–25 m, while `publication.py` requires `DesignBounds()` (59.5–2614 keV, 0.5–49.87 m).
- **Implementation:** reconcile the bounds; produce the NPZ by folding unit-fluence cells
  through `ArrivalResponseTable` (a folded unit cell is that cell's kernel; HPGe from a
  passing detector bake, NaI after its payload); a manifest builder around
  `env_scatter_export` + `ceelo_g4val --sphere-flux` (two independently seeded Livermore arms
  per anchor, 60 m world); run the registered-E0 gate.
- **Shows up:** `promote-candidate`.

### Environment GEANT4 evidence is integrated and near-field only
`[medium — validation]`
- **Effect:** the transport truth matches GEANT4 Livermore on the **integrated** scattered
  flux of five configurations (|z| ≤ 1.72), but not on the energy/angle shape a detector fold
  consumes; neither the open-air tier nor a concrete field is checked; nothing beyond 5 m
  although the campaign reaches 25 m (ŝ re-steepens past ~25 m); the "scoring sphere is a
  void" approximation is unmeasured. Published albedo/build-up data are sanity bands, not
  gates (Chilton–Huddleston is dose albedo and ANS-6.4.3 build-up is exposure-weighted —
  fold with μ_en/ρ or they read 1.5–2.5× off).
- **Implementation:** compare (E, μ) maps against `--sphere-flux` on the five existing
  configurations; add a 25 m field anchor (R = 50 cm on both sides), one air-tier and one
  concrete-field anchor, and one detector-present-vs-void run — via the manifest builder
  above.
- **Shows up:** quoting environment shape accuracy or distances beyond ~5 m.

### Rooms and clutter are not in the compact environment model
`[medium — feature gap]`
- **Effect:** a room can be computed node by node (direct room MC agrees with GEANT4 to
  −0.57%, z −1.72), but there is no reusable room model, and superposed single-surface
  kernels fail (DESIGN.md, Previously tried): 18–34% of room scatter comes from chains
  touching two or more surfaces. Vehicles, benches, adjacent containers and collimator walls
  are the same physics one level out. A single "clutter" scalar (GADRAS's effective wall
  distance) cannot represent the multi-bounce fraction and discards geometry the user often
  knows.
- **Implementation (designed and prototyped):** bake only a geometry-free per-material
  **slab-albedo operator**: φ_room = direct + air kernel + first-albedo quadrature over the
  room's FINITE faces + a patch radiosity/Neumann iteration for chains across surfaces. The
  prototype closes to −1.0 to −5.5% over five rooms (122 keV included; superposition
  +19.3%) and tracks a +48% near-wall swing to ~3%, at ~1 MB and ~25 min of bake per material,
  MC-free in <1 s per (room, E0); mixed materials, asymmetric rooms and off-centre detectors
  come free. Suggested home: `src/transport/SlabAlbedo.{h,cpp}` +
  `src/io/EnvRoomSolver.{h,cpp}`; non-convex rooms are out of scope. A one-number clutter
  fallback, if offered, should be a documented degradation of this model.
- **Shows up:** indoor measurements; near vehicles, large structures or collimators.

### Coincidence corrections under environmental scattering
`[medium — coverage]`
- **Effect:** summing factors use per-gamma detection probabilities in source + detector
  geometry only. Environmental scatter adds a large scattered continuum to the total
  efficiency that drives summing-out (0.45–1.76× the direct flux) while barely moving FEP, so
  a measurement over ground or in a room has larger summing-out than the vacuum answer CeeLo
  reports. Unlike attenuation, scatter is not separable once a scattered photon from one
  gamma can coincide with a full-energy other.
- **Implementation:** fold the environment kernel into the per-gamma total-efficiency term of
  the summing estimators, keeping FEP on the direct path; the two gammas share the geometry
  but are otherwise independent, so the total term likely factorizes while leaving a
  correlated scattered+full-energy residual — quantify it against a correlated-gamma GEANT4
  run (`--correlated-gamma`, air/ground world) first. Needs the composition layer.
- **Shows up:** cascade nuclides at standoff over ground or in a room.

### Importance-sampled environment primaries are not implemented
`[low — bake cost]`
- **Effect:** `PrimaryMode::VertexSampled` is declared but throws; primaries are analog.
  Histories to 0.3% at 662 keV grow steeply with distance (4M at 1.9 m, 26M at 6.9 m, 36M at
  13.4 m, >100M at 26 m), and a larger scoring sphere does not help the estimator.
- **Implementation:** sample the first vertex from an importance density on the prolate halo
  around the source–detector line (design in `NextEventEstimator.h`), accepted only through
  the paired-bias and variance × wall-time gates.
- **Shows up:** far-field nodes; if the adaptive rounds cost too much.

---

## Full-spectrum: source and charged-particle models

### Commit the Nov 3 source-escape and charged-particle payloads
`[medium — decision]`
- **Effect:** the Nov 3 2026 re-bake (`bake_work/fu/results/`, uncommitted; the previous
  payloads in `data/full_spectrum/` predate the shield x-ray components, and the old charged
  payload no longer loads) passes every source-escape and charged-particle gate. The
  detector payloads did not publish (TODO items above).
- **Implementation:** decide whether a gate-passing family is committed on its own or only
  with a complete bake, and whether `bake_work/fu/results/` replaces `data/full_spectrum/`.
  Before committing the charged payload, re-select the in-flight yield structure
  (`charged_particle_structure.json`, chosen on earlier banks) on the new banks, scoring the
  holdout planes and the off-grid final test together.
- **Shows up:** anyone wanting the fitted source or charged models.

### Full-spectrum payloads stop at Z = 92
`[low — coverage]`
- **Effect:** the photon data reach Z = 98, but the fitted payloads and the model's Z
  machinery stop at 92: `SourceModelPayload` and `ChargedParticleSourcePayload::
  max_effective_Z`, the Gaussian effective-Z blends over `Z = 1..92` in
  `FullSpectrumModel.h`, `PhotonuclearCrossSectionRow::mass_attenuation_cm2_g[92]`, and the
  campaign drivers' `Z > 92` validation. A Pu
  shield is outside every payload's validity range.
- **Implementation:** widen the blends and validation to `kMaxZ`, add Z 93–98 cells, widen the
  photonuclear row, and re-bake source and charged payloads together (the `Z/92` feature
  scale is part of a payload's meaning, and widening moves predictions near Z ≈ 90–92).
- **Shows up:** Pu/Am/Cm/Cf sources or shields.

### Distributed birth depth: activity in the shield, SNM self-attenuation
`[medium — coverage]`
- **Effect:** `SourceGeomFluxModel` assumes the primary enters the shield at depth 0 with full
  intensity (every base integrates a vertex over depth `u` against `exp(-tau u)`), so it
  cannot describe a source distributed through its own attenuating matrix — the dominant real
  case being U/Pu self-attenuation (the 186 keV U-235 line and the Pu 375/414 keV group are
  strongly self-absorbed, and their ratio to harder lines is how enrichment/mass is
  inferred), plus activated or contaminated shielding, NORM in concrete or soil, and waste
  drums. Two cases must both work: the emitter is the matrix (SNM metal/oxide), and a trace
  nuclide in an unrelated matrix — so the kernel must be parameterized on the matrix, not the
  nuclide — and either must hand off to ordinary outer shields. (The MC side is the "more
  source shapes and activity distributions" item.)
- **Implementation:** a double integral over (birth depth, vertex depth) changes the escape
  kernel, so it is a new component, not a new coordinate. Check first whether the
  uniform-birth case has closed forms the mechanistic bases can absorb (the uncollided term
  becomes `(1 - exp(-tau))/tau`), in which case the fitted residual may transfer.
- **Shows up:** safeguards / enrichment / holdup measurements; contaminated structures.

### Curved shields: the slab-geometry assumption
`[medium — accuracy]`
- **Effect:** the shield-scatter surrogate integrates vertices over a planar depth coordinate,
  i.e. a laterally infinite slab. Curved, finite shields (a sphere or cylinder around a
  source, a drum wall, a pipe) change the escape solid angle from a vertex and the scattered
  leg's chord distribution, so the same areal density scatters differently; today that error
  is absorbed into the fitted residual at the training geometry.
- **Implementation:** measure first — GEANT4 at fixed areal density in slab / spherical-shell
  / cylindrical-shell geometry (the spread *is* the effect). Then carry t/R (thickness over
  radius of curvature, signed convex/concave, derivable from the `SourceGeometry` shell stack)
  on the mechanistic bases, or failing that fit a correction factor. Not a port of GADRAS's
  scalar slab factor.
- **Shows up:** thick shields on small radii (t/R ~ 1): drums, SNM spheres, pipe walls.

### Source-secondary residuals
`[low — accuracy]`
- **Effect:** (a) deep-τ bremsstrahlung on the off-grid holdout is 4.5% median but 34.6% p95;
  (b) `pair_annihilation_continuum` is the worst-conditioned component (conditional-yield span
  8.4× on hydrogen-free rays, 12.4× with hydrogen); (c) `compton_electron_brems` is the
  weakest brems channel (span 4.4×, 7.2× at Z ≤ 20, 9.5× with hydrogen) because its parent
  energy is the *maximum* Compton recoil while the real recoil spectrum is broad;
  `pair_annihilation_in_flight` was accepted at slope +0.249 against the 0.25 pre-fit bar
  and read +0.260 in the next realization (now its own bound, 0.31). Held-out
  **element** yield error is 5.6% median (the 7-D coordinate thins the node cloud for element
  extrapolation). The acceptance grid's p95 tail (10.4%) is 59–250 keV opaque or high-Z cells
  with real errors (e.g. E59 Z82 1 g/cm² +21%, z 7.9).
- **Implementation:** (c) weight the log-k average by the actual recoil spectrum; (a)/(b) more
  bank statistics at τ ≥ 4; more trained elements or a per-axis bandwidth for Z
  extrapolation.
- **Shows up:** brems channels at τ ≳ 3, low Z; untrained elements; opaque low-energy shields.

### Compton-scattered shield x-rays are not modelled
`[low — spectrum fidelity]`
- **Effect:** the source-escape model emits the shield's K and L x-rays as atoms
  (unscattered); those that Compton-scatter on the way out form a continuum below each line
  that no component carries. On 73 bank cells (Z 50–92, 59 keV–12 MeV, τ 0.2–8, 2e5
  histories, ≥ 400 line counts) it is 2.2% of the lines at the median, p95 3.7%, max 4.9%
  (U, 3 MeV, τ 0.2), rising with τ (Pb at 100 keV: 0.3% at τ 0.2 to 4.5% at τ 8). Against
  everything escaping in the band it occupies (Kα backscattered to Kβ) it is 1.6% median,
  p95 3.3%, max 4.1%. In lead 80% of it lies within 10 keV below Kα, so at NaI resolution it
  reads as a ~2% line deficit, at HPGe resolution as a low-side shoulder. The off-grid score
  reports the share (`offgrid.fluorescence_continuum_share`). The charged-particle model
  likewise leaves scattered x-rays in its bremsstrahlung.
- **Implementation:** not planned while spectrum fits free the shield x-ray amplitude (the
  lines and their downscatter both depend on the geometry, and the free amplitude absorbs
  the NaI case). If the HPGe shoulder matters: a per-family continuum, the line's base times
  a fitted yield and a single-Compton downscatter shape of the family's lines.
- **Shows up:** the region just below the Pb K lines behind thick lead.

### Brems shape floor below ~120 keV in thick high-Z shields
`[low — accuracy]`
- **Effect:** the fitted brems shape puts yield where the truth is ~0 (E3000 Z82: 0–25 keV
  model/MC 12.9× at τ = 3, 40.9× at τ = 4) — 1.5–2% of the brems yield, ~0.1% at a caller,
  but most of that cell's +2.9% total over-prediction.
- **Implementation:** fit weighting (a per-channel floor or a log-space shape objective) and
  more bank statistics below 120 keV. Not interpolation, basis resolution or an envelope
  (DESIGN.md, Previously tried).
- **Shows up:** τ ≳ 3 high-Z shields, bottom ~4% of the escaping energy range.

---

## Full-spectrum: composition, runtime and tooling

### No payload loader; `full_spectrum_nuclides` ignores the source payload
`[low — usability]`
- **Effect:** `examples/full_spectrum_nuclides.cpp` loads only the detector and charged
  payloads; `payload.source` stays the default mechanistic model, so the only end-to-end
  histogram example does not show a fitted source-escape payload. Every caller assembles
  `FullSpectrumPayload` by hand.
- **Implementation:** a `--source-payload` option and a loader for the
  `data/full_spectrum/<family>/` layout (or one `CEELO_FULLSPEC` file per detector).
- **Shows up:** anyone judging the functional model by the example.

### `predict` folds the scattered source continuum on 16 uniform bins
`[medium — spectrum shape]`
- **Effect:** `predict` reduces the source measure's continuum to `source_quadrature_bins`
  (default 16) uniform bins over [0, E0], evaluates the detector once per bin at its centre,
  and spreads the tracking atoms over the bin as a flat box, so the scattered photons'
  full-energy deposits render as a staircase. Uniform bins also ignore where the structure is:
  at 662 keV the [165, 207] keV bin straddles the single-scatter back-scatter edge at 184 keV
  (the source measure itself is smooth at 1 keV resolution, with that edge and the
  forward-scatter shoulder below E0 its only sharp features). In the CeeLo-vs-GEANT4 study
  (10 cm, HPGe), 64 bins cut the shape error against the MC (TV on 4 keV bins) from 2.5 to 0.8%
  for 122 keV behind 5 cm PE, 2.2 to 1.6% for 1332 keV behind 2 cm Fe and 1.8 to 1.5% for
  662 keV behind 0.5 cm Fe; 256 add nothing but cost. The source evaluation is the expensive
  part, ~14 ms of CPU per source bin (9 s for 1 keV bins at 662 keV), against ~0.2 s for 16.
- **Implementation:** (1) evaluate the source continuum at landmark-aware nodes: breaks exactly
  at the back-scatter edge E0 / (1 + 2E0/mc²) and at E0, geometric refinement toward both,
  ~32–48 nodes; (2) treat the result as a density, monotone-cubic in each interval between
  breaks; (3) fold it without bins: the tracking atoms deposit density × the detector's yield
  per channel (yields interpolated between log-spaced incident nodes), and the detector
  continuum is a hat basis over those nodes, cached per (channel grid, cos θ), as
  `ArrivalResponseTable` already does for scattered fields. Then only physical edges remain.
- **Shows up:** shielded or scattering sources: back-scatter peaks and the scattered photons'
  full-energy continuum.

### `predict` is far slower than the 10 ms aspiration
`[low — performance]`
- **Effect:** 3.8–4.9 s CPU for 100 lines × 4096 channels (~400× the target); at 16384
  channels Na-22 takes 0.5 s, Am-241 4.5 s, Th-228 21.5 s, Eu-152 31.5 s. ~90% is the source
  model's per-channel cross-section and angular lookups (`interpolate_log_lin` on the angular
  curves alone ~45%), even at zero areal density; the detector response is ~7%. For fitting
  this is the binding cost.
- **Implementation:** an early out for a bare source; evaluate the smooth scattered continua
  on a coarse internal grid and rebin; cache cross sections and angular norms per line;
  lines in parallel; cache per-line source measures when only activity or resolution move.
- **Shows up:** fitting spectra with many lines or fine channel grids.

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

### IC-electron sum-peak magnitude is unresolvable against the current reference
`[low — validation]`
- **Effect:** `test_cascade_summing::ic_electron_g4_magnitude_ba133_2cm` asserts the sign of
  the ON-vs-OFF shift of the Ba-133 `sum356+81` peak (a sum peak gains what the singles peaks
  lose) and bounds the magnitude with the ON/G4 band, because which flag state is closer to
  GEANT4 (ON/G4 1.009, OFF/G4 0.998) is a coin flip against the reference's own ~1%
  statistics.
- **Implementation:** a ~25× larger Ba-133 GEANT4 ion run (~37M decays, air world,
  long-decay threshold) and ~36M MC decays per flag state — several hours to settle a ~1%
  question on a small peak.
- **Shows up:** quoting the IC-electron channel's accuracy on coincidence sum peaks.

### GEANT4 references are still scored at the old 1.5 keV FEP window
`[medium — validation]`
- **Effect:** `kDefaultFepWindowKeV` is 0.75 keV, but the GEANT4 references for configs 1–12
  and 25/26 were generated at 1.5 keV, as were the C++ suite's numeric expectations, so the
  gate and tests pin 1.5 keV (`kTestFepWindowKeV` in `tests/test_fep_window.h`, used by the
  suite and `cascade_ref_common.h`; `kReferenceFepWindowKeV` in
  `examples/benchmark_mc_configs.cpp`). Until the pins go, they do not measure the
  library's default. Configs 27/28 are already at 0.75 keV on both sides; the harness and
  `generate_all_spectra` default to 0.75 keV (pass the harness `--fep-window 1.5` to compare
  with the older references). The gap is real where near-peak forward scatter is: config 8 at 59.5 keV
  differs by −3.09% (z −13.5); bare NaI does not move.
- **Implementation:** settle the default (the MC has no resolution smearing, so the window
  is a choice of how much near-peak scatter to credit), regenerate the references at it,
  re-baseline the suite, and delete both pins and the note in `tests/data/geant4_reference/README.md`.
- **Shows up:** comparing CeeLo-at-default against the references; configs 7 and 8 most.

### GEANT4 physics-list mismatch between the references and the environment tier
`[open — validation]`
- **Effect:** every detector reference is `G4EmStandardPhysics_option4`, while the harness
  default and the environment tier are Livermore. The lists differ only in the low-energy
  Compton model, an angular (not attenuating) difference: 2.19% at 122 keV and 0.39% at
  2614 keV on a pure-scatter observable, 0.44% on the uncollided term. FEP- and
  attenuation-dominated results are insulated; total efficiency's scatter-in, continuum shape
  and source leakage channels are exposed, and 0.39% at 2614 keV is comparable to the pooled
  ≥2 MeV FEP residual (~+0.4%). (The EPICS2023-vs-2017 data axis is documented in DESIGN.md.)
- **Implementation:** run one point-source config under both lists (`--em-option4` vs the
  default) and compare FEP and total: ≲0.3% closes this as documented; ~1% means re-making
  the references under Livermore or carrying a stated systematic. Do it with the GEANT4
  K-label anchor regeneration.
- **Shows up:** scatter-dominated comparisons against the committed references.

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

### W-wall box total reads +2–3% (config 21)
`[low — residual]`
- **Effect:** config 21 (NaI behind a 10×15×20 cm W-wall box of cellulose) reads total
  +2.0% at 662 keV and +2.9% at 2614 keV (z 3.8 / 5.6) against its 16M option4 reference,
  inside the gate's tolerance + 3σ. It read +1.2% / +1.3% before the walls emitted K x-rays
  and source positrons were transported; both are required physics (the same fluorescence
  agrees with GEANT4 to ~1% behind Pb shells in the full-spectrum study) and move the Pb
  box (config 20) the right way (−0.96% → +0.26% at 662 keV). A/B: fluorescence adds ~+1%
  at 662 keV, the positron change +0.5–0.7% at 2614 keV.
- **Implementation:** a W-specific excess the missing fluorescence was masking. Check the W
  K-line set (EADL line energies and Kα/Kβ ratio against GEANT4's) and the escaped W K
  x-ray spectrum in both codes; regenerate the config-21 reference at 64M to halve its σ.
- **Shows up:** totals behind tungsten collimators and boxes.

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

### cfg 8 Marinelli FEP deficit at 59 keV
`[low — residual]`
- **Effect:** the water Marinelli runs FEP −1.9% (z −6.4) at 59 keV and −0.6 to −0.8% at
  100–200 keV on this branch's engine; totals ≤0.4%. The GEANT4-harness u/s decomposition localises it: the
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
  shell/density corrections) or a cfg-3 reference offset. In-crystal brems and Molière escape
  each move FEP 6–12% at 3 MeV; the single/double-escape ratio separates them (it responds to
  brems, barely to Molière; `--no-crystal-brems` / `--no-moliere` A/B them).
- **Shows up:** LaBr3 photopeaks ≥2 MeV, at the +0.5% level.

### Escape peaks and the bare-crystal 511 keV line at 2614.5 keV
`[low — residual]`
- **Effect:** with the pair cross section arbitrated (EPICS2023 = XCOM = GEANT4 Livermore;
  DESIGN.md "Known differences"), what is left at 2614.5 keV is not cross section: there
  EPICS2023, XCOM and Livermore agree to 0.3% and option4 sits ~1% above them for Ge–I (so
  a cross-section effect would make CeeLo read low). In the full-spectrum study the 1 cm CZT's double
  escape reads 1.023 ± 0.002 (z +14.8) against option4, HPGe single/double escape
  +0.9 ± 0.2% / +1.3 ± 0.3% (z +4.1 / +4.9) against Livermore (−0.4% / −1.0% against
  option4), and the HPGe net 511 keV line +4.7 ± 0.8% (z +6.2) against both lists.
- **Implementation:** separate the positron's annihilation point (in flight vs at rest,
  range) and the 511 photons' escape: tally the in-flight fraction and the annihilation-depth
  distribution in both codes for the bare CZT and HPGe at 2614.5 keV (see "In-crystal
  positron annihilation: residual approximations").
- **Shows up:** double-escape and 511 keV peaks of thin crystals at Tl-208 energies.

### HPGe backscatter peak at 59.5 keV reads 1.5–1.75× GEANT4
`[low — spectrum fidelity]`
- **Effect:** in the full-spectrum study (`scratch/fullspec_g4`, `hpge_60_bare`: Detective X
  fit crystal in the Micro-Detective housing) CeeLo / option4 is 1.48 at 46–47.5 keV and 1.75
  at 47.5–48.3 keV (z +19, +21), falling to 1.06 just above 48.28 keV, the 180° Compton
  energy of 59.54 keV; the K-escape lines beside it read 1.06. It predates the Nov 2026 fixes
  (1.7–1.9 before them) and the follow-ups did not move it. Only this scenario shows it:
  behind 1 cm polyethylene 1.05, HPGe at 122 keV 1.05, NaI 1.02 in the same window; the
  12–46 keV continuum of the same spectrum reads 1.10 (z +14). At 59.5 keV nothing
  crosses the crystal, so the peak is photons that pass beside it, backscatter off the
  housing behind and re-enter through the back or side.
- **Implementation:** split the CeeLo and GEANT4 runs by the volume of the last scatter
  (housing parts vs crystal) and compare the back/side passive layers of the CeeLo
  geometry with the exported GDML before suspecting the low-energy Compton model.
- **Shows up:** the 44–48 keV region of HPGe Am-241 spectra (the band is 1.2e-3 of the FEP
  in GEANT4; CeeLo's excess 5e-4).

### L x-rays below 10 keV and M x-rays are not emitted
`[low]`
- **Effect:** L lines are transported only at or above the 10 keV cut (Z ≥ 69); below it
  every L vacancy, and every M vacancy, deposits locally, because the photon tables start
  at 10 keV. On a bare 1 cm CZT, GEANT4 shows Cd/Te L-escape lines at E₀ − 3.1…4.3 keV
  (~0.02% of the FEP) that CeeLo lacks; with the Detailed electron walk they are what is
  left of the shelf just below the line: CZT's E₀−10 keV…E₀ region reads ~0.5 against
  GEANT4 at 60–662 keV (full-spectrum study).
- **Implementation:** extend the photon tables (at least photoelectric and total
  attenuation) below 10 keV, lower the L line cut in `generate_eadl_relaxation.py`, and
  extend the L1/L2/L3 photoelectric-share table to the lighter elements.
- **Shows up:** thin detectors of Cd/Te/I (CZT, CsI, NaI surfaces), low-energy spectra.

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
  ratios against EADL (the line energies are now EADL binding-energy differences, matching
  the measured lines).
- **Shows up:** just above high-Z K edges behind thick high-Z attenuators.

---

## Cascade / true-coincidence summing

### Na-22 511 keV drifted from GEANT4 after the positron-model change
`[medium — accuracy]`
- **Effect:** re-running `examples/cascade_observables` at the reference statistics gives
  **Na-22 P511 +3.24% (z +23.8)** and sum1785 +3.81% (z +6.0) against GEANT4 (everything else
  |z| ≤ 3.2; engine/G4 for P511 moved 1.0012 → 1.0340). It bisects to d3e0ee3, which replaced
  `emit_annihilation` with a full β⁺ spectrum + positron-range + in-flight-annihilation model;
  the GEANT4 reference did not change, so the drift is engine-side. It is still inside the
  gate's band and went unnoticed because the cascade block is not in ctest (Validation item
  above).
- **Implementation:** triage the positron model against GEANT4 for Na-22 first, then
  regenerate `our_cascade_multi.csv` (~85 s) and record the producing commit — regenerating
  first would bake the drift in as the expectation.
- **Shows up:** β⁺ emitters' 511 and sum peaks.

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

### Side & back incidence: what the side & back near field leaves open
`[low — accuracy in corners of the side & back regime]`
- **Effect:** the side & back near-field tables (stage 5) answer close beside and behind the
  detector, but:
  - in the shadow regime (direct transmission s < 0.3) the removal-only kernel misses full-energy
    in-scatter, which varies too steeply between nodes (20-70% off measured), so it stays gated
    (`Shadowed` / `NeedsMc`), table or not;
  - a source inside the assumed body (`back_body_depth_cm`, 15 cm in InterSpec), and every source
    beside or behind a far-field-profile response, is still carried by the kernel alone
    (`NearFieldUnmodeled`);
  - theta is measured from the crystal-face origin, so the back face of a long detector (a 2x4x16
    in NaI log) occupies the last few degrees, and its far-field back nodes there sit within ~2a of
    the housing;
  - boxes reuse the phi = 0 tables at every azimuth;
  - InterSpec's imported grids keep their behind-the-face rows only as far-field eta nodes, and
    their Monte-Carlo total donor runs without the side & back near field.
- **Implementation:** for the shadow regime, much denser nodes there or a scatter term in the
  kernel; for long detectors, place back nodes on the envelope's own surface rather than by angle
  from the face.  `closed_loop` probe families and transfer-mode back anchors could follow.
- **Shows up:** sources against the back of a thick housing at low energy; long logs; box crystals.

---

## Electron / CSDA physics

### Impact ionization is K-shell only, and low against option4 at 100 keV in CdZnTe / NaI
`[low — spectrum fidelity]`
- **Effect:** the Detailed and source walks ionize the K shell only (EEDL MF=23 MT=534); an
  electron makes no L vacancies, so Pb L x-rays from β sources in lead come only from
  photoabsorption. In electron beams with the cap at 100 keV (`CEELO_EWALK_EII_EMAX`, as
  option4) the E₀ − E_Kα line matches option4 at 300–1000 keV (Ge 1.04–1.10, CdZnTe
  0.95–1.00, NaI 0.98–0.99, |z| ≤ 2.1), but at 100 keV CdZnTe reads 0.77 ± 0.04 (z −5.2)
  and NaI 0.70 ± 0.04 (z −6.8), against Ge 1.05 ± 0.02. Option4 ionizes with PENELOPE's
  cross sections there.
- **Implementation:** compare EEDL MT=534 with PENELOPE's (Bote–Salvat) K-shell cross
  sections for Cd, Te and I at 30–200 keV, and run 100 keV GEANT4 `--pixe` beams (Livermore
  subshell tables, from an older EEDL; at 300 keV they agree with CeeLo within 1.02–1.08)
  to separate the cross section from the walk's depth distribution. For L, add
  the EEDL L1–L3 subshell tables (MT=535–537) and relax them with `sample_l_relaxation`.
- **Shows up:** E₀ − E_K lines of thin high-Z detectors; L x-rays of β sources in lead.

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
  above Z = 92 (EEDL2023 MF=26 MT=527 would be a new lineage), and the EEDL K-shell
  electroionization table above Z = 92 (today Np–Cf have no impact ionization). Cheaper
  for the radiative part: rescale uranium's by Z(Z+1)/A (error −2.5 .. +3.9%).
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

### Source-electron brems residuals
`[low — research]`
- **State:** primary Seltzer–Berger photon rates use the absolute normalization; condensed steps
  sample the Poisson photon multiplicity; energetic Møller/Bhabha delta electrons are transported for
  their bremsstrahlung. At 2.29 MeV over Z = {6, 26, 48, 92} and shield depths of {0.2, 0.5, 2}
  electron ranges, escaped yields are generally within ~2% of higher-statistics GEANT4 anchors.
- **Remaining:** the mean primary path is ~2–3.5% short in thin Fe/Cd spheres and ~3% long near the
  0.5-range Fe/Cd cells. This is the residual Highland geometry/path conversion, not the photon-rate
  normalization. The 50 keV delta-electron transport floor undercounts *created* secondary brems by
  10–30%, but escaped >10 keV photons from that component agree much better (absolute residual below
  ~0.001 photon per primary in the audited cells). A straggling-only A/B was rejected because it caused
  4–5% contained-U deficits.
- **Implementation:** validate at additional incident energies before adding any Z/depth correction;
  ultimately replace the Highland mean-displacement approximation with an electron-specific Mott/G-S
  step-to-geometry transform. Lower the delta floor only if HPGe-resolved comparisons show a visible
  low-energy residual.
- **Shows up when:** mono-electron/beta brems from thin elemental-equivalent source shields.

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

## Deferred by scoping decision

The working energy regime is **≤ 3 MeV**; items here are real but off the critical path.

### Photonuclear: parked; if revived, go per element
`[deferred]`
- **State:** giant-dipole absorption has its threshold at the least-bound nucleon separation
  energy (~6–8 MeV heavy, 16–19 MeV light), so below ~4 MeV it is 10⁻⁴–10⁻⁶ of the
  attenuation; it reaches 1–2% of escaping photons only near 12 MeV in high Z. The
  source-escape payload can carry a fitted term (`PHOTONUCLEAR=1`; evaluation off by default):
  0.59% median at trained elements, 13.4% off-grid, 27.2% held-out element.
- **Why per element:** leave-one-element-out over the ten trained elements is 23.1% median
  (p95 257%) — secondary multiplicity follows nuclear structure, not Z, so no effective-Z
  interpolant helps. The escaping secondaries also carry sharp nuclear de-excitation lines
  (Pb at τ = 0.2: 796.9 / 562.5 / 890.6 / 2601.6 / 1757.8 keV at +27 to +45σ) that must be
  atoms from isotope-resolved data, not a more flexible continuum (≤1.9% of escaping photons,
  <0.5% of the spectrum). Until then GEANT4 is the truth for prompt nuclear lines.
- **Cost:** the `PHOTONUCLEAR=1` GEANT4 bank is ~450 CPU-hours (1.1–4.0 ms per event at
  8–12 MeV); a cheaper recipe needs a fresh accuracy study first.

## Minor

- Rename the `MCDET_NO_{STRAGGLE,B2,TAILS}` env knobs (`ElectronCsda.cpp`) to `CEELO_`.
- `PhotonTransport.h` still mentions `diag_interactions`, now a dev-only study tool.
- `ElectronCsda.h`'s `kSB_min_energy_keV` comment says the brems tables start at 200 keV;
  they start at 10 keV.
- Doppler broadening omits the F(p_z) profile-reweighting factor (few-% tail) and the
  relaxation cascade from the Compton vacancy — Compton-edge **shape** only; never cite it
  for FEP/total.
- Full-spectrum: K-line yields are gated off above their support limits (320 / 175 / 280 /
  150 keV), leaving a ~0.02%-of-FEP step — a decay envelope outside the node support would
  remove it.
- Full-spectrum: the secondary surrogate's leave-one-interior-τ-node-out bandwidth CV doubles
  the local gap and so over-widens (on-grid yield error 3.24% vs 1.19% on-node); a fold that
  brackets the query like the real off-grid cells would fix it — judge on the acceptance grid
  and off-grid holdout.
- Environment `run-streaming` is packer-bound (18.6 s/node against 8.0 s of transport):
  overlap the next node's transport with packing.
- Environment `fit-archive` crashes on reduced designs when a component has no supported
  validation node (it stays in `model_inventory`): drop skipped components from the
  inventory or fail with a message naming the component.
- `env_kernel_bake --landmark N --out FILE` text kernels drop the per-node energy edges
  (`--binary-out` is correct): write/read them in text too, and reject a mismatched `SHAT`.
- `surface_flux_sigma` (face / side / back) still sums per-bin variances in quadrature (×1.485
  low for the total); `ChainClass` labels are noisy in room corner wedges.
- Fit caches skip `require_full_verification` in `publication.py`, and the env README says
  `promote-candidate` needs a full seal while `workflow.py` accepts any valid seal: pick the
  rule and make code and README agree.
- `refine_multiple_scatter_cells.py` re-runs a refined cell at its CSV-rounded energy, rewrites
  unchanged model files and is not idempotent; `fit_multiple_scatter_surrogate.py` silently
  skips a cell whose model file is missing.
- Detector-response scoring differs run to run in the last digit (≤1.6e-14 relative,
  summation order); payloads are byte-identical.
