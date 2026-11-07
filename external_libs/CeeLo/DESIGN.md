# CeeLo — design and validation notes

Developer-facing record of how CeeLo is built, what physics it models, how it is
validated against GEANT4, and what is deliberately approximated. For an
introduction to the library see [README.md](README.md); open problems live only in
[TODO.md](TODO.md) (each entry records its **Effect**, what its **Implementation**
would take, and where it **Shows up**; entries are deleted when resolved). This
document records what *is* modelled and why.

## Purpose

CeeLo is a C++17 Monte Carlo photon transport library for computing gamma-ray and
X-ray detector efficiency. Given a detector geometry (NaI, LaBr3, HPGe, CZT, ...),
optional attenuators/shielding, and a source configuration, it computes:

- **Full-energy peak (FEP) efficiency**: probability that a photon deposits 100% of its energy in the crystal
- **Total efficiency**: probability that a photon deposits any energy in the crystal
- **Pulse-height spectrum**: energy-deposit histogram (compatible with spectral analysis software)
- **Cascade summing corrections**: coincidence summing-in and summing-out factors
- **Full-spectrum approximation**: fitted, MC-free operators that reproduce those spectra in
  milliseconds for arbitrary source, shield, distance and environment (below)

The goal is agreement with the large, fully detailed Monte Carlo engines (in practice
GEANT4 11.4.0) to ~1% from 30 keV to 3 MeV, with anything beyond that understood and
documented (Validated Configurations and Known differences below). The full-energy
peak, the total efficiency and the pulse-height spectrum are all primary outputs: the
full-spectrum approximation is fitted to CeeLo's labelled spectra, so continuum shape,
escape peaks and x-ray lines have to be right too, not just the two efficiencies.

## Build notes

Build commands, dependencies and the CMake options are in [README.md](README.md#build).
Use `-DCMAKE_BUILD_TYPE=Debug` only for assert checking: it is several times slower,
and the MC test suite (~15–30 min in the default `RelWithDebInfo`) grows with it.

**Portability** (both learned getting the vendored copy to build with GCC 12 on Linux —
neither shows up on macOS/libc++, so do not "clean them up" without a Linux build):

- The library links `Threads::Threads` **PUBLIC**. `EfficiencyCalculator.cpp` and
  `ResponseGenerator.cpp` use `std::thread`, which on libstdc++ needs `-pthread` at link
  time; on macOS pthread lives in libSystem, so the omission is invisible there.
- RNGs are seeded through the `seeded_rng()` helper in `EfficiencyCalculator.cpp`, which
  uses `std::seed_seq`'s **iterator** constructor. libstdc++ rejects the `initializer_list`
  form for `uint64_t` arguments. Its three words are chosen to reproduce the streams the
  `initializer_list` form produced, bit for bit — adding or reordering words would move
  every seeded result in the library, including the GEANT4-gated validation numbers.

**rapidxml** (only for `CEELO_BUILD_RESPONSE`) is NOT vendored — for ODR safety the host
application's own copy is used. Configure searches a few sibling `3rdparty` directories
for `rapidxml/rapidxml.hpp` (SandiaDecay and SpecUtils each bundle a copy); otherwise pass
`-DCEELO_RAPIDXML_DIR=<dir containing rapidxml/>`. The cascade adapter's SandiaDecay
source is located the same way (`CEELO_SANDIADECAY_DIR`); point it at the host's copy
when vendoring.

## Repository layout

The tracked file set is the reusable library: `git archive HEAD` produces exactly the
payload that gets vendored into a host application. CeeLo's own efficiency results are
regenerated rather than committed (`tests/data/ceelo_reference/README.md`), paired 1:1
with the GEANT4 references in `tests/data/geant4_reference/`;
`profiling/compare_validation.py` gates them against each other.

```
src/
  cascade/          True-coincidence summing: analytic + correlated-MC estimators,
                    beta/IC spectra, optional SandiaDecay adapter
  cross_sections/   Compiled-in photon data (generated element_data.cpp etc.) and lookup
  efficiency/       EfficiencyCalculator - main entry point
  export/           GDML and macro export for GEANT4 validation
  full_spectrum/    Full-spectrum approximation: source-escape, charged-particle and
                    detector-response operators and their composition
  geometry/         Ray tracing (cylinder, box, sphere), attenuator shells, source geometry,
                    environment geometry (air + ground/room slabs); GeometryDescriptor,
                    the storable detector geometry with a dependency-free text codec
                    (its XML codec is in io/)
  io/               Storable detector-response payload, MC-free response kernel,
                    EFFTRAN-style efficiency transfer (CEELO_BUILD_RESPONSE);
                    environment-scatter kernel tables, transport blocks and fitted payload
  materials/        Material definitions (composition, density), builtin materials
  physics/          Electron CSDA range tables, stopping powers, Moliere scattering, FEP window,
                    2BS bremsstrahlung angular density (BremsAngular)
  transport/        Photon transport (Compton, PE, Rayleigh, PP, fluorescence, bremsstrahlung);
                    next-event estimator, FluenceTally, Rayleigh angular evaluator
examples/           Examples; benchmark_mc_configs + bench_configs.h (validation configs)
tests/              Boost.Test suites + reference data
tools/
  prepare_cross_sections/  Locked EPQ/EPICS2023/ESTAR generators and reports
  geant4_validation/       GEANT4 harness (ceelo_g4val); generate_all_spectra (built with examples)
  full_spectrum/           Full-spectrum truth campaigns, fitters, gates and bake drivers
                           (README.md is the bake runbook; env_scatter_fit/ is the
                           environment-scatter campaign package)
  efftran_transfer_demo, make_golden_response, response_presets.h
profiling/          compare_validation.py — the GEANT4 gate
visualization/      Dependency-free WebGL viewer for exported GDML geometries
```

One-off studies, parameter sweeps, external-code benchmarks and design notes live in a
`studies/` tree that is **not** part of the distributed source; the root `CMakeLists.txt`
adds it only `if(EXISTS studies/CMakeLists.txt)`. A few source comments name those
studies for provenance; their conclusions are summarised here.

## Key Design Decisions

### Physics

- **Transport**: every photon interaction (Compton, PE, Rayleigh, pair production) is
  sampled explicitly; histories are analog unless a biasing technique applies (see Event
  Biasing — forced first interaction is auto-enabled for optically thin crystals).
  Electron energy deposition uses CSDA with Molière multiple scattering for KE > 200 keV.
  The crystal's condensed-history walk applies the Highland log term on the accumulated
  path (not per sub-step), adds Bohr energy-loss straggling to each step while the
  electron's KE exceeds 500 keV (the final, budgeted step stays deterministic), and walks
  until KE → 0 rather than dropping a residual when the step budget runs out. Those three
  changes (Aug 2026) closed a coherent −0.4 to −0.5% FEP deficit at ≥2 MeV and most of the
  thin-crystal (CZT) electron-escape deficit; a matched GEANT4 process-isolation
  decomposition exonerated the Seltzer–Berger brems sampling (≤0.2%). Env knobs
  `CEELO_NO_CRYSTAL_STRAGGLE` / `CEELO_DIAG_STRAGGLE_ALL_KE` A/B the straggling.
  This is the **Fast** walk (`ElectronTransport::Fast`, the default): straight below
  200 keV and without large-angle scattering above, so near a surface it neither
  diffuses nor backscatters. The **Detailed** walk (spectrum fidelity, below) fixes the
  partial-deposit continuum that costs; FEP and total do not need it.
  A walk that leaves all geometry reports `ElectronDepositResult::escaped`: an escaping
  positron emits **no** annihilation photons (the world outside is vacuum; CeeLo does not
  transport charged particles in air), where it used to annihilate at the last surface
  it crossed.

- **Cone importance sampling**: for point sources with no shielding, photons are emitted
  within a cone subtending the detector, each event weighted by the cone's solid-angle
  fraction (~100× speedup for typical geometries). The estimator is
  `(1/N) * sum(w_i * I_i)` — the denominator is N (events), NOT the sum of weights.

- **Cross-section data**: photon magnitudes and angular factors come directly from locked
  **EPICS2023 EPDL** records for **Z = 1–98 (H–Cf), 10 keV–20 MeV** (`kMaxZ`,
  `kPhotonDataMinEnergy_keV`/`kPhotonDataMaxEnergy_keV` in `CrossSectionData.h`; outside
  that window the accessors return the endpoint value, so a caller that must not clamp
  should refuse the query). Each process has its own per-element grid with source nodes,
  thresholds and both sides of absorption edges; adaptive probes gate the maximum material
  representation error at 0.1% and integrated error at 0.05%, and angular tables also
  gate CDF and mean-cosine errors. Generated tables must not be hand-edited; regenerate
  and audit them through `tools/prepare_cross_sections/`. Against NIST XCOM,
  photoelectric + Compton + pair agree to ≤0.3% for Z 92–98 at 15 keV–10 MeV and ≤0.5%
  for every element compared at 12–20 MeV. The coherent term differs: for Z 92–98 EPDL's
  is 1–34% below XCOM's, mostly just below K edges (the anomalous-scattering correction
  EPDL applies and XCOM omits), moving the total *including* coherent by up to −2.6% (U)
  and −4.8% (Am, 122 keV). Electron-side tables stop at uranium (see Known differences).

- **K-shell fluorescence**: after photoelectric absorption, K vs outer shell is sampled
  from EPICS2023 EPDL MT=534/MT=522 fractions; K vacancies produce fluorescence X-rays or
  Auger electrons from EPICS2023 EADL MF=28/MT=533 relaxation data through Z = 99
  (radioactive-decay daughters one beyond the photon tables). Line energies are EADL's
  binding-energy differences E_b(K) − E_b(filling subshell), not its tabulated transition
  energies: against xraylib's measured lines the intensity-weighted error is 2.4 eV over all K
  lines (tabulated: 131 eV; Pb Kα1 74.970 keV against 74.969 measured, where the tabulated
  value was 75.250), and each line records its filling subshell (`line_filling_subshell`,
  which sets the Kα/Kβ family and the L vacancy left behind). The same K relaxation
  (`sample_k_relaxation`) follows photoelectric absorption in **source material and
  shields** (`set_source_fluorescence`, default on; the x-ray joins the source secondaries
  as `Origin::Fluorescence`) and, with spectrum fidelity, the K vacancy a Doppler-sampled
  **Compton** collision leaves (`TransportConfig::enable_compton_fluorescence` in the
  detector; in source material and shields whenever their Compton collisions are
  Doppler-sampled, which spectrum fidelity and the source-escape campaign turn on).

- **L-shell fluorescence** (lines ≥ 10 keV, so Z = 69–99: Pb, W, Th, U, ...): an L vacancy
  — from photoabsorption outside the K shell (L1/L2/L3 by the EPDL MT=535–537 shares of the
  photoelectric cross section, tabulated per element with duplicated nodes at every edge),
  from the subshell a K-alpha's electron left, or from a Doppler-sampled Compton collision —
  relaxes through `sample_l_relaxation`: a retained line with the subshell's yield, else
  Coster–Kronig transfer (f13, f12, f23) and another try, else Auger. The photoelectron
  carries E − E_Li. Lines below the 10 keV cut and every L vacancy of a lighter element stay
  local (the photon tables start at 10 keV), so HPGe, NaI and CZT histories are unchanged;
  CZT's Cd/Te L-escape lines are not modelled (TODO.md). The cascade vacancy model takes a
  K-alpha's L vacancy from the same filling subshell.

- **Bremsstrahlung**: Seltzer–Berger scaled differential cross sections (Z = 1–92, the
  complete NIST EPQ set; Np–Cf reuse uranium's), **27 native electron energies from
  10 keV–20 MeV × 32 κ**, `float32`, from pinned public-domain NIST EPQ. κ = 0.025 is a
  quadratic through EPQ {0, 0.05, 0.1} and κ = 0.075 a cubic through {0, 0.05, 0.1, 0.15};
  everything else is copied from equal-κ EPQ nodes, with no row normalization or
  element-specific correction. Rejection sampling uses a 1/k envelope; brems photons are
  transported recursively. Aggregate GEANT4 spectra/rates are kept only as validation
  observables.

- **Compton scattering**: Butcher–Messel sampling (as GEANT4's `G4LowEPComptonModel`)
  with S(x,Z) bound-electron rejection. (An earlier Kahn's-method implementation had a
  subtle bug worth +5–7% FEP.) The scattered energy is Doppler-broadened by the bound
  electron's momentum (impulse approximation, analytic subshell profiles).

### Spectrum fidelity (`enable_spectrum_fidelity`)

`EfficiencyCalculator::enable_spectrum_fidelity(true)` switches on the detector physics
that shapes the pulse-height spectrum but moves FEP and total by at most a few tenths of a
percent, so FEP/total runs leave it off (the default) and full-spectrum work turns it on
(the detector campaign's `--detector-spectrum-fidelity`, which the bake passes;
`fullspec_ceelo` in `scratch/fullspec_g4`). It sets three `TransportConfig` fields and two
`SourceGeometry` switches, each with its own setter for A/B runs:

- **`enable_compton_fluorescence`** (`enable_compton_fluorescence`): the K vacancy a
  Doppler-sampled Compton collision leaves (`ComptonScatterResult::shell_index`) relaxes
  through `sample_k_relaxation` and its x-ray is transported as after photoelectric
  absorption (`FluorescenceRecord::from_compton`); off, the whole binding energy is
  deposited on the spot. It restores the K-escape peaks above ~300 keV (NaI, CZT: 0.90–0.95
  of GEANT4 → 1.00–1.02).
- **`electron_transport = ElectronTransport::Detailed`** (`set_electron_transport`): every
  crystal, attenuator and source electron goes through
  `ElectronCsda::deposited_in_scoring_detailed` instead of the Fast walk — a mixed
  (class II) condensed-history scheme in the PENELOPE mould at every energy:
  - *Elastic:* per material, a Wentzel distribution whose screening and normalization
    reproduce the first two transport cross sections of EPICS2023 EEDL MT=525
    (`electron_elastic_data.cpp`, 92 elements × 36 energies, 26.6 KB; across EEDL's
    256 keV–10 MeV tabulation gap each moment is interpolated as its ratio to screened
    Rutherford). Collisions above the cut are sampled one by one with mean free path
    λ_h = max(λ_W, min(C1 λ₁, C2 R)), C1 = C2 = 0.1; the softer ones become one Gaussian
    deflection at a uniformly sampled point of the step (random hinge).
  - *Inelastic:* collisions transferring more than W_cc = 2 keV are sampled one by one from
    the Møller cross section; each delta ray is walked like any electron from where it was
    born, and the primary takes the binary-collision recoil. The soft remainder (ESTAR CSDA
    loss less the hard collisions' mean) is drawn from PENELOPE's artificial distribution
    with the soft variance. Brems (Seltzer–Berger, `sample_step_bremsstrahlung`) and
    positron in-flight annihilation are as in the Fast walk. Positrons use the electron
    elastic moments and Møller collisions.
  - *Boundaries and cost:* a flight that reaches a region boundary stops there and a new
    step is sampled beyond it. An electron whose CSDA range is below
    `Geometry::safety_distance` (a lower bound on the distance to every surface of the
    crystal, active volume, bore and attenuators) deposits in place without walking, and a
    flight shorter than the safety distance skips the ray trace.
  - Env knobs for A/B: `CEELO_EWALK_C1`, `CEELO_EWALK_C2`, `CEELO_EWALK_WCC` (0 = CSDA loss
    with Bohr straggling above 500 keV and no delta rays), `CEELO_EWALK_EABS` (absorption
    energy, 1 keV), `CEELO_EWALK_NO_FASTPATH`.
- **`enable_impact_ionization`** (`enable_impact_ionization`, which also sets
  `SourceGeometry::set_source_impact_ionization`): electron-impact K-shell ionization from
  the EPICS2023 EEDL K-subshell electroionization cross section (MF=23 MT=534,
  `electron_ionization_data.cpp`: EEDL's native nodes from threshold to 20 MeV, linear-linear
  as EEDL declares, Z ≤ 92, 15 KB) times ω_K. In the Detailed walk each flight draws a
  Poisson number of K x-rays with mean rate × path (the CSDA yield of the rest of the
  slowing-down for a flight that stops, and for a contained electron in place; tabulated per
  element on a log grid and zero at and below its threshold, which falls between grid nodes),
  the element by its share and the line by intensity; the x-ray leaves isotropically from a uniform
  point of the flight as a brems photon does (`BremsPhoton::impact_Z`, recorded as
  `FluorescenceRecord::from_impact`). The stopping power already holds the ionization's
  energy, so the x-ray is paid out of the walk's scoring deposit (`unpaid_impact_keV` when the
  electron left first). Source and shield electrons emit them per step of their walk, and an
  electron that is not walked emits its slowing-down yield at its birth point
  (`Origin::ImpactIonization`). `CEELO_EWALK_EII_EMAX` caps the electron energy (100 keV
  mimics GEANT4 option4 without PIXE, whose PenIoni ionizes only below 100 keV).
- **Source and shield Compton Doppler** (`SourceGeometry::set_compton_doppler_enabled`): their
  Compton collisions sample the bound electron, so the K and L vacancies they leave relax
  like the detector's.

Electron pencil beams into 2 cm slabs against GEANT4 option4 (1 µm cut;
`scratch/fullspec_g4/ebeam`), full-energy fraction CeeLo/GEANT4, Fast → Detailed:

| | 30 keV | 50 keV | 100 keV | 300 keV | 1 MeV |
|---|---|---|---|---|---|
| Ge | 1.50 → 1.015 | 1.50 → 1.019 | 1.50 → 1.022 | 1.21 → 0.998 | 1.13 → 0.998 |
| CdZnTe | 1.70 → 1.009 | 1.72 → 1.010 | 1.73 → 1.013 | 1.32 → 1.006 | 1.22 → 1.000 |
| NaI | 1.63 → 1.010 | 1.66 → 1.013 | 1.67 → 1.012 | 1.28 → 0.992 | 1.20 → 0.993 |

The Fast walk has no backscatter below 200 keV at all; the Detailed walk's partial-deposit
bands (fractions of E₀) sit within ±7% of GEANT4. Impact ionization restores the K-escape
line E₀ − E_Kα: capped at 100 keV as option4 is, CeeLo/GEANT4 for it reads Ge 1.05 ± 0.02 /
1.04 ± 0.03 / 1.10 ± 0.05 at 100 / 300 / 1000 keV (without it 0.10 / 0.40 / 0.77) and CZT
(Cd) 0.99 ± 0.03 at 300 keV (0.79); the 90–99% band of a 100 keV Ge beam 0.32 → 0.90. The
full-energy fraction at 30–100 keV still reads +1% (Ge 100 keV +1.5% → +1.0%), so impact
ionization was only part of that excess. Uncapped (the default) the line is 4.2× option4 at
Ge 300 keV and 2.0–2.4× in CdZnTe and NaI at 300–1000 keV, where option4 has no impact
ionization; GEANT4 with PIXE on (its Livermore subshell tables, from an older EEDL) agrees
with it at 300 keV: Ge 1.04 ± 0.01 (z +3.1), CdZnTe 1.08 ± 0.02 / 0.97 ± 0.03 (Cd / Te line),
NaI 1.02 ± 0.02, and Ge 1.07 ± 0.02 at 1 MeV, with every deposit band within 5% (the 100 keV
cap reads 0.16–0.51 there).
At 100 keV, where both have it, CdZnTe and NaI read 0.77 ± 0.04 and 0.70 ± 0.04
(z −5.2, −6.8) against Ge's 1.05 ± 0.02: option4 ionizes with PENELOPE's own cross sections
below 100 keV, not EEDL's (TODO.md). Step size is converged: C1 = C2 = 0.01, 0.05 and 0.1
agree within ~2% band by band; 0.2 starts to drift. W_cc = 1 and 2 keV agree. Without the Møller
collisions (CSDA loss) the walk diffused correctly but its backscattered electrons were too
soft and, at 1 MeV, too few (early backscatter −20%, full-energy fraction +2%); in photon
spectra that put a 1 cm CZT's FEP at 2614.5 keV 1.8% above GEANT4 (now −0.5%).

In photon spectra (`scratch/fullspec_g4`, bare crystals) the fidelity run matches GEANT4
Livermore: HPGe FEP within ±0.1% at 662–2614.5 keV, Compton continuum, valley and
backscatter within ±0.25% at 356–2614.5 keV; with impact ionization (32M events) the 1 cm
CZT reads +0.15 ± 0.06% / −0.08 ± 0.08% at 356 / 662 keV (z +2.6 / −1.1) and 1332 keV
+0.03 ± 0.11% (before it), NaI 3×3 at 662 keV +0.02 ± 0.04% (z +0.4), HPGe 662 keV
−0.01 ± 0.06%. Against option4 FEP reads +0.1–0.6% at 356–1332 keV: option4 itself sits
0.4–0.9% below Livermore in FEP there, a GEANT4-list difference. Impact ionization leaves
FEP, total and the 10–100 keV band unchanged within their noise and raises the crystal
K-escape peak E₀ − Kα against option4 (no PIXE): NaI 0.98–1.01 → 1.07–1.08 and CdZnTe
0.96–1.01 → 1.04–1.10 at 356–2614.5 keV, HPGe 1.00–1.18 → 1.21–1.74 (z +2.4 to +6.6), the
part option4 does not emit; behind 0.2 cm Pb at 2614.5 keV the Pb K lines read 1.2× option4
(shield electrons' impact ionization).

Cost (CPU time, single thread, same binary, fidelity on/off): HPGe 1.3× at 59.5–662 keV
and 2.6× at 2614.5 keV; NaI 3×3 0.7× at 662 keV (the fast path skips the Fast walk's ray
traces in a large crystal) and 3.2× at 2614.5 keV; 1 cm CZT 2.0× at 662 keV and 6.3× at
2614.5 keV. The default (Fast) path is unchanged.

### Distances, frames and reference points

Every distance in CeeLo is one of three things, and confusing them is the easiest way to
introduce a silent few-percent bias that no internal self-consistency test catches.

- **The crystal frame.** The one frame the tracer works in. `z = 0` is the front face of
  the **crystal solid**, `z = L` its back face; a source in front is at **negative z**.
  Radial distance is from the symmetry axis.
- **The dead layer is carved out of the inside of the crystal.** The solid still spans
  `[0, L]` at radius `R`; only the *active* volume shrinks, to `[t_front, L - t_back]` and
  radius `R - t_side`. So the crystal face stays at `z = 0` however thick the dead layer
  is, and **a dead layer must never be added to an outward distance**. "The front of the
  crystal" and "the front of the active volume" are planes `t_front` apart — say which.
- **Only attenuator shells sit in front of the crystal.** Each extends its front thickness
  toward negative z, so the "endcap front" (the face a user can touch) is at
  `z = -Σ(layer front thicknesses)`.

`ReferencePoint` selects which plane a caller's `dist_cm` is measured from, and
`GeometryDescriptor::endcap_front_offset_cm()` is the **only** conversion between them,
applied in `DetectorResponse::query_position()`. Anything else building a source position
from a user distance must go through it. (It once added the dead layer, placing every
`EndcapFront` source ~0.7 mm too far on a typical HPGe — invisible far-field, a ~3%
deficit at contact.)

Two nearby functions look like the same bug and are **not** — do not "fix" them:

| function | includes dead layer? | why |
|---|---|---|
| `endcap_front_offset_cm()` | **no** | a real distance, query-time only, so correcting it costs nothing |
| `transverse_half_extent()` | yes, deliberately | a *scale parameter* computed identically at generation and query time; changing it desynchronises stored responses |
| `Geometry::outer_bounding_radius()` | yes, deliberately | a sampling *bound*; over-estimating is safe, under-estimating loses paths |

Rule of thumb: a quantity used at generation **and** query is a model parameter and must
not change; one used only at query time is an interpretation and can be corrected.

### Geometry

- **Detector shapes**: cylinder and box, declared with `set_detector(material,
  CylinderDims{radius, full_length})` or `set_detector(material, BoxDims{half_x, half_y,
  full_length})` — transverse extents are HALVES, the axial extent is the FULL length.
  That mix is deliberate and is what `GeometryDescriptor::dimensions_cm` serializes; the
  convention is stated once, in the CRYSTAL DIMENSION CONVENTION block of
  `src/geometry/Geometry.h`. (Sources use half-extents throughout: `half_length`,
  `half_dims`.)

- **Bulletized (rounded) front edge**: `set_bullet_radius(r_b)` replaces the cylinder's
  outer front corner with a quarter-torus fillet (ANGLE's `bulletizingRadius`); 0 is a
  sharp edge. The ring circle has radius `rho_c = R - r_b` in the plane `z = r_b`.
  Removed volume (Pappus): `V = 2 pi rho_c r_b^2 (1 - pi/4) + pi r_b^3 / 3`. The solid is
  convex and a subset of the sharp cylinder, so `Cylinder.cpp` computes the sharp interval
  and moves only an endpoint that lands inside the removed corner onto the fillet:
  safeguarded Newton on `G(t) = (rho - rho_c)^2 + (z - z_c)^2 - r_b^2`, convex on every
  corner-zone sub-segment; grazing chords that stay in the corner zone minimise the same
  `G`. The closed-form ray–torus quartic is deliberately not used (catastrophic
  cancellation at near-double roots — exactly the grazing case); it survives as a comment
  and a test oracle. Every consumer (photon transport, the electron walk, `ResponseKernel`,
  the stochastic-Rayleigh re-trace) goes through `Geometry::trace_ray`, so this covers all
  of them; cone biasing uses the bounding cylinder, which bulletization only shrinks
  inside. Attenuator shells stay sharp (the endcap does not follow the fillet). With a dead
  layer the active volume gets the inward-offset fillet, radius
  `r_b - max(t_front, t_side)`. Sharp crystals pay one predicate and run the original code
  bit for bit.

- **Rounded bore tip**: `set_bore_hole(r, depth, rounded_tip=true)` caps the bore with a
  hemisphere of the bore radius (ANGLE's `<core rounded="yes"/>`), keeping the stated
  depth (the apex sits where the flat bottom would). The capped bore is convex within its
  own z range, so its ray interval is the hull of the tube's and the ball's — closed form.
  For a blunt tip (`depth < 2 r`) the ball reaches past `bore_z_end`, so
  `rounded_bore_interval` clips it to `z <= bore_z_end`.

- **The dead layer is `solid - active - bore`**: a bore is empty over its whole length,
  including where it runs through dead material (a back dead layer, a bore reaching the
  front dead layer, a bore wider than the active radius, or the removed corner of a
  bulletized crystal). `trace_cylinder_geometry` cuts the bore from the dead pieces as well
  as from the active volume (at most 3 dead pieces). When the whole bore lies inside the
  active volume the cut is skipped and the original code runs bit for bit. InterSpec's
  `DetectorGeometryDiagram` (`v_dead = v_solid - v_bore - v_active`) and `Geant4Export`
  (one polycone whose `rmin` is the bore) describe the same solid.

- **Stored geometry (`GeometryDescriptor`)**: the serializable form of all of the above
  and the only thing a saved response carries; `build_geometry()` and
  `ResponseGenerator::configure_calculator()` replay it. `bullet_radius_cm` and the bore's
  rounded-tip flag are written to XML only when non-default (`bulletRadius`,
  `roundedTip`), so sharp, flat-bored responses keep their exact bytes and
  `content_hash`. `sm_xmlSerializationVersion` is deliberately **not** bumped for purely
  additive attributes — the only version logic is a range check, so a higher version
  would be rejected outright by existing builds. `Geometry`'s preconditions are
  `assert`s; `GeometryDescriptor::problems()` is the release-safe mirror (plus too-short
  `dimensions_cm`, a dead layer that consumes the crystal, a bore wider than the *active*
  radius, a bore clearing the outer fillet but not the active one). `build_geometry()` and
  `configure_calculator()` run it and throw; callers building geometry from user input
  should consult it first and relax rather than throw (InterSpec's ANGLE `.outx` import
  drops an unrepresentable fillet or rounded tip with a warning).

- **Attenuators** (detector shielding): cup-shaped shells around the crystal,
  `add_attenuator(material, front_thickness, side_thickness, z_start, z_end)`, concentric
  (each wraps the previous) and always traversed.

- **Back structure.** A layer whose `z_end` reaches past the envelope it wraps is a cup closed
  by a back plate of that excess thickness - no extra field, so older builds trace it the same.
  `GeometryDescriptor::set_layer_back_thickness_cm` / `layer_back_thickness_cm` work in that
  convention (setting one moves every outer layer's back, and a collimator that reached the old
  back, so their own plates are unchanged); `back_structure_extent_cm` is how far the described
  structure reaches behind the crystal (vacuum gaps included; 0 = nothing described), and
  `outer_envelope()` the detector's true outside. A vacuum layer with only a back is a gap
  behind the crystal. The GDML export subtracts the inner envelope whenever a layer has a front
  or a back.

- **Source shielding**: spherical shells (point source) or concentric shells matching the
  source shape — `add_source_shield(material, thickness)`, `(material, t_radial, t_end)`
  (cylinders), `(material, t_x, t_y, t_z)` (boxes). Per-dimension thicknesses apply to both
  ± faces; components may be zero (not all). Call `set_*_source()` first. Full mode
  transports in the shields (Compton/PE/Rayleigh with re-tracing); FEP-only mode uses
  stochastic Rayleigh.

  Full-mode source transport also generates **secondaries** for every source shape
  (`SourceFullTransportResult::secondaries`): K x-rays from photoelectric absorption in
  source material/shields; 511 keV pairs from pair production there, emitted where the
  positron stops (or in flight), and none for a positron that escapes the source
  (`set_source_positron_transport`, default on; with source-electron transport on, an
  escaping positron is handed to the detector, where it annihilates when it stops);
  bremsstrahlung from electrons stopping there (Seltzer–Berger along a
  Molière walk, `ElectronCsda::walk_in_source_geometry`); and escaped electrons that reach
  the crystal (`enable_source_electron_transport`). The Molière diffusion matters:
  straight-line CSDA overestimates escape from dense shields 2.4–10×. These channels
  matter above the 1022 keV threshold and for recoil electrons >200 keV; for config 11
  (0.5 cm Fe) at 3 MeV they add ~+4.2% (annihilation), ~+3% (brems) and ~+1% (electron
  escape) of total efficiency. `set_source_brems(false)` and
  `benchmark_mc_configs --no-source-brems` / `--no-source-electrons` A/B them. Electrons
  whose CSDA range cannot reach their birth region's surface
  (`SourceGeometry::min_distance_to_boundary`) skip the walk and use the geometry-free
  `sample_brems_in_material` (same physics, escape provably impossible; ~1.2× faster on
  configs 11/12).

- **GDML export**: attenuators export as cup-shaped subtraction solids (outer full
  cylinder minus a shorter inner cylinder) so the front face is covered (a hollow tube
  once cost −26% at 30 keV). The export mirrors `RayTrace`'s accumulation exactly: the
  first shell sits on the crystal radius (the dead layer is carved *inside* the crystal, so
  it adds nothing there), and a shell gets its subtraction solid whenever it reaches past
  the cavity it encloses — front *or back*, so a layer whose `z_end` runs past the back face
  exports its back plate. A shell of zero side thickness is just its front/back disks and is
  written as plain tubes (a union for both), avoiding a subtraction of two equal-radius
  cylinders. A detector with a dead layer exports the whole crystal as a non-scoring
  `CrystalLV` carrying the active region (`ActiveSolid`, `active_crystal`) as a daughter in
  the same local frame — the tracer's dead = crystal − active − bore, with the active fillet
  radius `r_b − max(front, side)` — so no booleans are involved. Every GDML also carries a
  `<!-- ceelo_cone_bounds_cm r= zmin= zmax= -->` comment: the outer bounding cylinder CeeLo
  sizes its own emission cone from, which the harness's `--cone-bias` uses. A crystal with a **bore and/or bulletized edge** is one
  native `<polycone>` (`CrystalSolid`): `rmax` follows the fillet arc, `rmin` is the bore
  radius (via the tip hemisphere if rounded); arcs become 64 chords (sagitta ≪ 1 µm). An
  earlier boolean-solid version produced ~3250 stuck-track warnings per million under cone
  bias (the polycone ~25, like a plain tube). The two profiles must be merged on the
  *union* of their z values, and the harness's `--cone-bias` bounds parser understands
  `<polycone>`. **The bore cavity must not get a daughter volume** (it lies outside its
  mother and G4 navigation sticks); it is left to the world material, which matches the MC
  when exporting with `vacuum_world=true`. GDML resolves `volumeref` only against volumes
  already read, so a daughter must be defined before the volume that places it.

### The full-energy window

`src/physics/FepWindow.h` holds the single definition of "in the peak":
`ceelo::kDefaultFepWindowKeV` (0.75 keV half-width, suited to HPGe). Everything that has
to agree reads it: the MC's FEP tally, the FEP-only early kill
(`TransportConfig::fep_window_keV`), the cascade `PeakWindow` default, and the analytic
in-window Compton fraction (`kn_in_window_fraction`).
`EfficiencyCalculator::set_fep_window_keV()` moves tally and kill together; the kill must
never assume a narrower window than the one scored, or FEP comes out low with no error.
`ResponseGenerator` records the window in `ResponseProvenance::fep_window_keV`. Callers
that know the peak override it: `CascadeSummingCalc` takes each window from the fitted
peak's sigma.

**Most committed GEANT4 references were generated at 1.5 keV** (configs 27/28 at
0.75 keV). Until the rest are regenerated, everything compared against them pins 1.5 keV:
`kTestFepWindowKeV` (`tests/test_fep_window.h`, the C++ suite and `cascade_ref_common.h`)
and `kReferenceFepWindowKeV` (`examples/benchmark_mc_configs.cpp`, which takes
`--fep-window` to override). The difference is real where near-peak forward scatter is:
config 8 at 59.5 keV is FEP 7.479e-02 ± 0.161% at 1.5 keV vs 7.248e-02 ± 0.163% at
0.75 keV (−3.09%, z −13.5); bare-NaI config 1 shows none (|z| ≤ 1.0). The harness and
`generate_all_spectra` default to 0.75 keV, so **pass `--fep-window 1.5` to the harness**
when comparing against a 1.5 keV reference.

**Rayleigh in the analytic FEP transmission.** `fep_survival_removal_mu` drops `mu_rs`
(coherent scatter cannot leave the window by energy), which treats a Rayleigh-scattered
photon as continuing undeflected — right for a thin layer, wrong for a thick one (in Fe at
60 keV half the coherent scatters exceed 20°, `<cos> = 0.80`). 
`rayleigh_deflection_loss_fraction(E, tau_nr, mat)` (`ResponseKernel.h`) gives the fraction
h(T) of deflected photons then absorbed, from the transport's own form-factor sampler
(h → P(θ ≥ 90°) as T → 0; 0.25 for 0.5 cm Fe at 60 keV). An analytic consumer adds
`mu_rs * h(mu_rem * t)` per layer OUTSIDE the emitting volume. This was all of InterSpec's
+10% model-over-MC residual behind 0.5 cm Fe at 60 keV. Approximations: h at the layer's
normal depth, per layer; the flat in-window Compton credit carries no deflection penalty
(`fep_survival_removal_mu_depth` is the empirical stand-in).

### Simulation Modes

- **Full mode** (default): tracks all photon interactions through attenuators and crystal,
  electron CSDA with bremsstrahlung; gives FEP, total efficiency and the spectrum.
- **FEP-only mode**: full MC inside the scoring crystal, but importance-weighted stochastic
  Rayleigh through non-scoring segments (attenuators, source shielding), treating
  Compton/PE there as absorption; photons leaving the crystal are killed. Faster; FEP only.

### Analytic cascade summing (`compute_cascade_analytic`)

`cascade/AnalyticCascade.h` gives per-window summing factors from a nuclide's
`DecayCascade`s and an injected efficiency provider, with no MC (the exact level-scheme
joint-probability form; see the header). Everything independent of efficiencies and
windows — level-scheme DAGs, x-ray line tables, summing-in occurrence lists — is prepared
once in an `AnalyticCascadeModel`; callers evaluating the same cascades repeatedly (a fit)
build it once. Each primary's summing-out factor and each summing-in pair key's coincident
survival are memoized, so a U-235/U-238 evaluation over ~300 windows takes tens of ms.
Providers must be deterministic and side-effect free.

The MC estimators (`compute_cascade`) can also deposit **internal-conversion electrons**
and the coincident K-Auger (`CascadeConfig::enable_ic_electrons`, default OFF, off path
bit-identical): the conversion electron is transported through the source geometry and
air and deposited only when it can physically reach the crystal (an automatic
distance/solid-angle/air gate). It matters at contact/well/4π (Ba-133's 356 keV summing
factor ~0.44 → 0.31 at contact) and matches GEANT4 per-decay peak areas to <1% at 2 cm.

### Decay x-ray models (SandiaDecay adapter)

`cascade_adapter::build_cascades` represents a branch's decay x-rays in one of two ways.

- **Vacancy model** (default, `CascadeOptions::vacancy_xray_model`): the branch carries
  its K/L vacancy sources (EC shell fractions, per-transition K/L ICCs) and the engine
  emits lines from the daughter's EPICS relaxation data. Line exclusivity per vacancy and
  per-decay multiplicity are exact, and with a valid level scheme the analytic engine
  enumerates x-ray–x-ray pairs exactly. The cost is ~25 line occurrences per converting
  transition.
- **Transition-group ("flat") model** (`vacancy_xray_model=false`, or no vacancy data):
  SandiaDecay's tabulated x-ray lines become `Xray` members linked to the branch's gammas
  at the gamma's intensity (independence). InterSpec's reference lines and sum-peak check
  use it for speed. With `kxray_coincidences` (default on) the K x-ray members are also
  linked to each other: with P_EC the K-capture probability, N_IC the K conversions and
  N_ICIC the pairs of K-converted transitions in one cascade, the expected vacancy pairs
  per decay are `D2 = P_EC*N_IC + N_ICIC`. N_ICIC walks the level scheme where both
  transitions are placed (each level splits its flux by occurrence weight); otherwise it
  falls back on SandiaDecay's gamma-gamma links (which often disagree with the scheme —
  Lu-177's are ~0.1×). A vacancy is filled by line a with
  `p_a = I_a / max(P_EC + N_IC, sum I_K)`; the pair (a,b) is emitted `D2*2*p_a*p_b` times
  per decay (stored as links P(b|a) = joint/I_a) and a line with itself `D2*p_a^2` times
  (`CascadeMember::self_pair_prob`). The analytic engine's summing-out of a flat x-ray uses
  these links plus `2*self_pair_prob/I_a`. Against the vacancy model's exact enumeration,
  D2 matches to machine precision where it is all IC×IC (Lu-177, Eu-155); where EC×IC
  enters the independence assumption gives +0.2% (Pb-203) to +5.9% (Te-123m).
  `build_cascades` cost: Eu-152 0.06 ms, the U-238 chain 0.7 ms.

### Aperture quadratures and stored-response identity (`CEELO_BUILD_RESPONSE`)

- **Quadrature storage.** An `ApertureQuadrature` keeps every ray's material stretches in one
  contiguous `segments` array (a ray holds an offset and count; iterate with `q.segs(ray)`), and
  each stretch carries a slot into the quadrature's distinct `materials`, so a kernel evaluates mu
  once per material per call rather than searching per segment. Build rays one at a time with
  `add_segment` (`discard_segments` drops a ray that turns out unused); rays traced in parallel into
  separate quadratures are joined with `append_rays`. Against per-ray vectors this made tracing
  1.3-1.7x and the per-energy kernel ~1.2x faster, bit-identical (Release, 2026-11).
- **`content_hash()`** identifies a response by exactly what `to_xml_string()` stores - the
  certificate and the version stamp excluded - but is hashed from the members (FNV-1a over raw
  values, length-prefixed, presence-flagged) rather than from the formatted XML, which cost ~1 ms a
  call. It mirrors `serialize_xml` field by field and condition by condition, so a new stored field
  must be added to both - and an optional one (written only when present or non-default) must
  contribute nothing to the hash when absent and a uniquely named `tag()` plus its content when
  present, so existing responses keep their hash, as they keep their bytes.
  `test_detector_response_io/content_hash_tracks_every_serialized_field` edits the first, middle and
  last token of every attribute and element value, and removes each element, one at a time, and
  fails if the hash does not change exactly when the stored content does - for the fields its test
  responses populate. The definition changed in 2026-11 (from hashing the XML text), so every
  response's hash changed once.

### Detector-side line sets / etendue (`CEELO_BUILD_RESPONSE`)

`src/io/DetectorEtendue.h` builds a fixed set of LINES through the active crystal, for
hosts that integrate an EXTENDED source. Reversing the order of integration of the
per-element kernel (`eps_fep_element`, Eq. 5) gives an integral over lines parameterised
on the DETECTOR side by a hull point `x` and direction `w`, with measure `dA |w·n| dΩ`
(`dV dΩ = dA |w·n| dΩ ds`). Everything detector-side of a line is independent of source
and energy, so one line set serves every energy and fit iteration, and the source
integral along a line is analytic for a uniform source. `build_etendue_lines` samples a
4-D Halton set (`src/io/LowDiscrepancy.h`) aimed at a target sphere;
`sample_hull_points` + `append_etendue_line` let a host use its own direction proposal.
`line_interaction_probabilities` (live cross sections) and
`DetectorResponse::fep_line_probabilities` (stored mu tables) give the per-line kernel.
Lines that miss are dropped WITHOUT renormalising (they count in `n`), which keeps the
estimator unbiased. Validated in `tests/test_detector_etendue.cpp` against the per-point
aperture kernel (~0.4% at 2^16 lines) and the analytic hull measure. InterSpec's
volumetric Activity/Shielding fit is the first consumer.

### Total efficiency close in (`CEELO_BUILD_RESPONSE`)

A generated `DetectorResponse` stores the total as a far-field tier (kernel, b(E) or
eta_tot(E,θ) at `far_distance_a`) carried inward by the ray-traced kernel. On a 5×5 cm HPGe
coax that holds to ≤1.3% on axis but drifts to 4.4% off axis close in (75°, 60 keV), so
the near-field nodes, which run full transport anyway, also keep their totals:
`TotEffPayload::near_field` holds ln(MC total / far-field model) on the FEP table's lattice
and applies it below its breakpoint (worst 0.44–0.94%, RMS ≤0.4% at 1–10 cm). A response
without it serializes exactly as before. A host keeping only the total and on-axis FEP
sets `GenerationOptions::scans_stop_on_total`, so the angular and near-field nodes stop on
total precision (6–16× fewer histories at ≥662 keV; worst total 0.84%).

### Side & back incidence (`CEELO_BUILD_RESPONSE`)

A source beside or behind the detector (theta > 90 deg from the crystal-face origin).
- **Generation:** `GenerationOptions::back_coverage` {None (default), Coarse 6 angles, Fine 11,
  93-180 deg} adds a far-field scan (stage 4, its own node ids, so front seeds are untouched) at
  `back_precision_scale` x the node precision, sampled with the general cone
  (`EfficiencyCalculator::enable_general_cone_sampling`: the outer envelope's bounding cone about
  its own axis, rather than the isotropic fallback whenever the source's z lies within the
  detector's). Its columns are prepended to `eta_fep` / `eta_tot` (cos theta < 0). FullMc
  fixed-grid only; `closed_loop` refuses it.
- **Side & back near field** (`GenerationOptions::back_near_field`, on by default; General and
  Contact profiles only): stage 5, its own node ids, measures `DetectorResponse::back_near_field`
  and its total twin `TotEffPayload::back_near_field` exactly as stage 3 measures the front table -
  ln N = ln(MC / (eta K)) at each node, an ln N = 0 anchor row - on (shape energy x the back-coverage
  angles x clearance), General 0.05 / 0.3 / 1 / 2 a, Contact eight rows 0.05-3.3 a.  The direction
  of the crystal's back rim is added as an angle when none lies within 4 degrees of it: close beside
  the detector the near field changes sharply as a source passes the crystal's back end (+8% at the
  rim of a 6.5 x 5 cm HPGe, between the 135 and 155 degree nodes, before it was added).  The distance
  axis is d = a + the CLEARANCE past `DetectorResponse::back_near_body()` (nearest point), not the
  distance from the face-centre origin, so no node can fall inside the detector and beside a long
  detector the axis says how near a source is; nodes are placed with
  `Envelope::position_at_clearance`.  The anchor row sits at the nearest far-field back node (where
  the back eta was measured, so ln N is 0 there); rows within 1.3x of it are dropped (long
  detectors).
  - **Assumed body** (`GenerationOptions::back_body_depth_cm`, recorded as
    `provenance.back_body_depth_cm`): the envelope's back face is moved to at least that z, so no
    node is spent where a real detector keeps its cryostat, PMT or electronics (InterSpec: 15 cm
    behind the crystal face).  It is not a placement limit: a source inside it - outside the
    described detector - is left to the kernel (`NearFieldUnmodeled`), as on a response without the
    tables, so there is a seam at the body's back face.  Dropped (0) when a far back node of the
    detector itself falls inside it.
  - **Shadow skip:** a node whose direct transmission is below `model_sigma::shadow_refuse_s` -
    where the evaluator refuses whatever the table says, and the MC would run to the CPU cap - runs
    no MC.  It, and any node with no full-energy counts or > 50% uncertainty, takes the nearest
    measured clearance at its energy and angle - else no correction - with a 100% sigma.  Never
    another angle's value: a shadowed column filled from beside the detector put a step at the
    body's back face that the line path's prefactor grid read as 0.5% (60 keV, GEM35-70).
  - phi = 0 only, as the front table; boxes reuse it at every azimuth.
- **Evaluation:** `EtaTable::eval_ln` interpolates the front and back node sets separately, so a
  covered response's front values, uncertainties and flags equal an uncovered one's, bit for bit -
  below the grazing node too, where both clamp to it.  Behind the face plane the gap from the last
  back node is bridged linearly up to the plane, where that clamped grazing value holds. Behind the face plane (`common_eval`): covered -> the usual floors; uncovered and
  not angle-flat -> `OutOfRangeClamped` + `model_sigma::behind_extrapolation`; angle-flat
  (curve / flat transfers) -> no clamp, the transfer envelope saturates (s2 = 1); behind the
  crystal's back plane with a described back -> + `back_structure_approx`, with none -> `NeedsMc`
  + `behind_plane` as before. A total table without back nodes under a covered FEP table is
  charged `behind_extrapolation` itself. The near-field gates use the distance to the nearest point
  of the outer envelope (`d_near_cm`), not the exit along the ray from the face centre, which
  overstates it beside a long detector. The collimator shadow gate (s = transmitted / geometric solid angle of
  the active crystal) applies in front for a collimated response, and behind the face because a
  removal-only kernel misses full-energy in-scatter once the direct path is mostly blocked. Against
  MC around two HPGe and two NaI detectors (approximate backs, 40-1332 keV):
  - without back nodes, behind the crystal's back plane the model read ~-11% at s 0.1-0.3, ~-33% at
    s 1e-4..1e-2 and ~-100% below, within ~5% above 0.3; beside the crystal, through thin side
    walls, within ~7% down to s 0.01 - so the gate applies behind the back plane (or a collimator);
  - with back nodes, within ~5-8% above s 0.3 but 20-70% off below it, beside the crystal too (the
    nodes' scatter-dominated values interpolate badly) - so the gate applies past the face plane.
  A back node the MC could not measure (no full-energy counts, or > 50% uncertainty) is held from
  its nearest measured neighbour with a 100% sigma.
  With that, `behind_extrapolation` = 0.05 leaves every region's RMS pull <= 0.77 over 2332
  points.
  With the side & back near-field tables, a query past the face plane and outside the assumed body
  applies the table (`ln_boost` at d = a + clearance past the body, its node sigma, no
  `NearFieldUnmodeled`) up to its anchor row; a query inside the body takes exactly the path above.
  Against the same truth (InterSpec `test_BackIncidence`, 2026-11-05: Coarse, General profile, 0.5%
  nodes, 15 cm body; s >= 0.3, outside the body), over GEM35-70, Detective-X and a 3x3 NaI: median
  |model/MC - 1| 0.5-1.0% and p95 1.5-2.4% behind and beside at clearances >= 0.3 a (without the
  tables 0.7-1.2%, p95 2.5-4.4%), and at contact (0.05 a) median 0.4-0.7%, p95 <= 3.1% (without:
  2.6-5.0%, p95 to 10%), RMS pulls 0.57-0.87 with no new sigma term.  The stage costs 4-10% of such a run (shadow-
  skipped nodes 4-14%).  A 2x4x16 in NaI log improves at the median (contact 4.2% -> 1.5%) but
  keeps a 17% p95 behind (see TODO.md). `back_structure_approx` (0.10) is the cost of a reasonable back description only:
  leaving 2 cm of steel out of an HPGe's back changed its efficiency directly behind 30-fold at
  122 keV.
- **Serialization:** a response with back nodes is stamped version 4, which older readers refuse
  (they would interpolate across the seam); everything else keeps its version and hash.  The side &
  back near-field tables are optional `<BackNearField>` elements (under the root and under
  `<TotalEfficiency>`) plus the `backBodyDepthCm` provenance attribute, written and hashed only when
  present; `finalize()` refuses one on a response without back eta columns.  A version-4 reader that
  predates them ignores the elements and answers as without the tables.
- **Clearance:** `geometry/ConvexSeparation.h` measures the distance between convex solids
  (`convex_separation`: a point, sphere, cylinder about any axis, or rotated box, by alternating
  projections) and where a ray enters one (`convex_entry_distance`).
  `GeometryDescriptor::Envelope::solid()` is the detector's outside as such a solid; place it with
  `ConvexSolid::transformed()`. A host uses these to keep a source and its shielding clear of the
  detector wherever it sits, and to find where a ray from the source reaches the housing.
- **Recorded choices:** `generate()` records `ResponseProvenance::back_coverage` (fixed grid),
  `precision_profile` (unless the caller passed its own `node_precision` map) and
  `n_anchor_angles` (transfer mode), so a host can offer the same run again. They are optional
  `Provenance` attributes, written and hashed only when set: a response saved before them loads
  with them empty and keeps its hash.

### Efficiency transfer (EFFTRAN-style, `CEELO_BUILD_RESPONSE`)

`src/io/EfficiencyTransfer.h` anchors an efficiency curve at ONE reference position (from
CeeLo MC **or** measured `(E, eff)` points) and transfers it via the ratio of the
interaction-weighted effective solid angle,
`eps(E, target) = eps_ref(E) · K(E, target) / K(E, ref)`, with
`K = ApertureQuadrature::interaction_omega` (the MC-free `ResponseKernel`). It is EFFTRAN's
T-ratio transfer, fully 3-D, Rayleigh-aware and chord-exact. Entry points: the standalone
`EfficiencyTransfer` class (`eps_at`, `eps_at_with_tsrc` for MEFFTRAN self-absorption) and
`make_transfer_response(...)`, which emits a normal `DetectorResponse`.
`GenerationOptions::transfer_mode` (+ `n_anchor_angles`) builds a DRF from MC at only the
on-axis backbone (~28 MC nodes vs ~544 for the General profile, 11–16× faster); a
measured curve needs no MC. `DetectorResponse::model_transfer` (`SigmaTransferModel`)
inflates the reported sigma off-axis and near-field (`NearFieldUnmodeled` below ~2a).
The total kernel credits passive-layer Compton scatter-in (`passive_compton_recapture`,
`kTotalScatterInRecapture = 0.8`), which took the HPGe contact total from −2.5% to −1.4%
with FEP unchanged. `tools/efftran_transfer_demo` reproduces EFFTRAN's `et.bat`. Tests:
`tests/test_efficiency_transfer.cpp`.

**Measured accuracy** (NaI 3×3 and an HPGe coax, far-field on-axis anchor vs direct MC):
FEP on-axis ≤1% for d ≥ 5a (median 0.2–0.55%), ~1.5–2% at 3–5a, ~3–3.6% median near
contact; off-axis ~3.5% median from one anchor, <1% median with `n_anchor_angles=3`. Total
on-axis ≤0.6%; the HPGe total still drifts off axis (use the few-anchor mode). Bottom line
for a host: a viable "quick DRF" and measured-curve → other-geometry estimator, ≤1%
on-axis far, a few % near/off-axis with honest inflated σ, full MC as the upgrade.

## Cone Biasing Policy

**Pure cone biasing is used only for point sources with no source shielding or source
material — or in FEP-only mode.**

With source shielding or an extended (non-zero density) source in full mode, photons can
scatter in the source and change direction, so pure cone emission would miss them. Those
runs use either isotropic 4π emission or, by default when the detector subtends
Ω/4π < 0.15 from the source centre (not Marinelli), the **two-stream estimator** (Event
Biasing), whose direct stream is cone-emitted and weighted by the analytic no-interaction
probability while its scatter stream stays isotropic. **Mixture angular biasing**
(`BiasingConfig::mixture_cone_alpha`) is also unbiased with source scatter but is only a
manual knob.

In FEP-only mode cone biasing is safe even with source material: the stochastic Rayleigh
weighting handles source effects without changing direction. Extended sources size the
cone per event from the sampled position.

A point source whose z lies within the detector's z range (beside it) gets no cone by default
- the cone is aimed along the detector axis, and from there it would be the whole sphere - so it
is emitted isotropically. `EfficiencyCalculator::enable_general_cone_sampling(true)` instead
bounds the outer envelope by a cone about the direction to it (`geometry/BoundingCone.h`), ~10x
fewer histories beside the detector; off by default, so fixed-seed results are unchanged.

The same rules apply to the GEANT4 harness (`--cone-bias` for point sources only;
`generate_all_spectra.cpp`'s `run_config(..., use_cone_bias)` sets both sides). Isotropic
GEANT4 runs need ~50–100× more events for comparable statistics.

## Event Biasing (Variance Reduction)

`BiasingConfig` (`EfficiencyCalculator::set_biasing()`) adds unbiased variance reduction.
When `set_biasing()` is NOT called, `compute()` enables techniques per geometry, energy and
termination targets (`compute_effective_biasing()`); a default-constructed `BiasingConfig`
forces fully analog sampling (used by A/B consistency tests).

- **Forced first interaction** (`force_detector_interaction`): the primary's first
  interaction is forced inside the detector; the history carries weight `1 − exp(−τ)`, the
  point drawn from the truncated exponential across segments. Secondaries are never
  re-forced. Off for the Marinelli primary path. For source-effect events the gate is per
  event: forcing applies only when the primary is the sole depositing branch. Wins where
  pass-through is likely (CZT ≥200 keV: FEP FOM ×1.2–1.4, total ×1.5–2.7); hurts FEP FOM
  in thick crystals. **Auto**: central-ray τ < 0.7 when FEP precision matters, τ < 4 when
  only total is targeted.

- **Mixture angular biasing** (`mixture_cone_alpha` = α): emit from
  q(Ω) = α/4π + (1−α)/Ω_cone with weight (1/4π)/q ≤ 1/α; isotropic draws landing in the
  cone get the in-cone weight; degenerate cones fall back to isotropic. FEP FOM ×5.6–29 on
  configs 11/12 but total FOM degrades, so it is **manual only**, superseded by two-stream.

- **Two-stream stratified estimator** (`two_stream`, `two_stream_direct_fraction` = f):
  for full-mode source-effect runs (non-Marinelli), every tally partitions exactly by
  whether the primary interacted in the source geometry
  (`SourceFullTransportResult::interacted`). A deterministic round-robin fraction f of
  events runs the **direct stream** — cone emission weighted by (Ω/4π)·T,
  T = `SourceGeometry::no_interaction_probability()`, then plain detector transport (forcing
  always valid). The rest run the **scatter stream** — isotropic analog source transport,
  killed when the primary did not interact. Weights 1/f and 1/(1−f) keep Σw·I/N unbiased;
  the reported variance is slightly conservative. The scattered stream still carries real
  FEP (small-angle Compton and Rayleigh-only chains stay inside the FEP window — up to
  ~24% of FEP for cellulose at 59 keV), so never assume FEP_s = 0. Measured with Compton
  bias γ = 0.3: joint min(FOM_fep, FOM_tot) ×2.9–6.1 (cfg 11), ×5.4–15.7 (cfg 12). f ≥ 0.85
  is dangerous at low energy (FEP_s undersampled). **Auto**: f = 0.5 when only FEP
  precision is targeted, else 0.25, only when Ω/4π < 0.15 (close geometries stay analog).

- **Compton-angle mixture biasing** (`compton_cone_fraction` = γ): at the primary's FIRST
  Compton vertex in the source geometry, sample the direction from
  q(ω) = (1−γ)·p(ω) + γ·1[ω ∈ cone]/Ω_d (p = KN × S(x,Z); the cone subtends the detector's
  bounding sphere), weight p/q ≤ 1/(1−γ); kinematics at the sampled angle are exact
  (`finish_compton_at_angle`). The angular norm is composite GL32 in s = sin(θ/2)
  (~1e-6, cached). Total FOM ×1.8–3.8 on shielded configs; γ ≥ 0.5 hurts large unshielded
  sources. **Auto**: γ = 0.3 unless only FEP precision is targeted.

Weight-composition rules (apply to any new technique):
- One weight per event, applied once at the tally; every recursive `transport_photon()`
  belongs to the same history.
- The estimator denominator stays N (`sum_weights += 1.0` per event).
- Branch-level Russian roulette is **unsafe for FEP**: the FEP indicator is non-linear in
  per-branch deposits, so killing a branch with survivor reweighting biases FEP even though
  it preserves the mean deposit. Only whole-event roulette/splitting is valid.
- FOM harness: `build/examples/benchmark_mc_configs --fom-csv <file>` with `--analog` /
  `--force-interaction` / `--mixture-alpha <a>` / `--two-stream <f>` / `--compton-bias <g>`
  / `--precision <p>` (FOM = 1/(σ_rel² · T)); every run also writes `our_<cfg>_multi.csv`
  (consumed by `profiling/compare_validation.py`). Consistency tests:
  `tests/test_biasing_consistency.cpp`. Stress configs 13–17, 19 in
  `examples/bench_configs.h` guard the auto-policy outside the GEANT4-anchored set.
- On source-electron configs the electron brems machinery is 40–60% of CPU; further FOM
  there is perf work, not estimator work.

## Source Geometry vs Detector Geometry

- **Detector geometry** (`Geometry` / `RayTrace`): the crystal plus attenuator shells. Ray
  tracing returns `PathSegment`s with material, thickness and an `is_scoring` flag.
- **Source geometry** (`SourceGeometry`): the radioactive volume plus its cores and
  shields (`set_source_material()`, `add_source_core()`, `add_source_shield()`), separate
  from the detector with open space in between.

The photon's journey: emission point → source material → source shields → air gap →
detector attenuators → crystal. Full mode uses analog transport in the source geometry
(`transport_source_photon()`); FEP-only mode uses stochastic Rayleigh weighting
(`compute_transmission_fep_only()`).

**Air between the source geometry and the detector is not transported.** An opt-in
analytic correction, `set_air_attenuation(AirAttenuation::AnalyticNoScatter)`, multiplies
each event's weight by exp(−μ·d) over the source-to-detector air path (μ without
Rayleigh); it is not yet validated against GEANT4. Air (and ground/room) scatter is the
separate, offline environment tier below.

### Source shapes (`SourceType` / `SourceGeometry::Shape`)

`Point`, `Cylindrical`, `Rectangular`, `Sphere`, `Marinelli`. Every extended shape
self-attenuates (`set_source_material`), takes concentric shape-matched shields, and
supports uniform or exponential-depth activity. Shape dispatch is centralised
(`sample_source_position`, `source_material_path`, `trace_source_segments`,
`compute_transmission[_fep_only]`, `min_distance_to_boundary`, GDML export), so `compute`
and `compute_cascade` pick up every shape for free.

- **Spherical** (`set_spherical_source(center, radius, rotation=I, inner_radius=0)`): solid
  ball or, with `inner_radius > 0`, a shell with a non-attenuating void centre.
- **Hollow / annular cylinder** (`set_cylindrical_source(..., inner_radius>0,
  inner_half_length=-1)`): a negative `inner_half_length` is a through-bore (closed-form
  annular sampler); a value shorter than the half-length is a closed cavity (rejection
  sampling) — what InterSpec's nested cylinders need. `inner_radius = 0` reduces
  bit-identically to the solid cylinder.
- **Hollow box shell** (`set_rectangular_source(..., inner_half_dims)`): outer box minus a
  centred inner void (rejection sampler; `box_shell_path` subtracts the void chord). See
  `tests/test_hollow_box_source.cpp`; pipe conventions in `tests/test_pipe_source.cpp`.
- **CylinderEndOn vs CylinderSideOn** (InterSpec orientations) are the
  `set_cylindrical_source` rotation: EndOn = identity; SideOn = 90° about y,
  `R = [[0,0,1],[0,1,0],[-1,0,0]]`.
- **Source cores** (`add_source_core`): attenuating layers filling **inward** from a hollow
  source's inner surface, given outermost first — the mirror of `add_source_shield()`.
  Cores never emit and never change the per-emitted-photon normalization; whatever is left
  at the centre stays a void. This is what lets CeeLo express InterSpec's
  Activity/Shielding model, whose source shell need not be the innermost layer.
- **One stored layer stack.** Cores, the emitting shell and shields are one
  **innermost-first** `std::vector<SourceLayer>` (`layers()`); "core" and "shield" name only
  which side of `source_layer_index()` a layer sits on. The ray march walks that vector and
  **the GDML export emits the same vector**. An attenuating interior makes segment ORDER
  load-bearing (`transport_source_photon_impl` places interactions by accumulated length),
  so those sources take an ordered concentric march; with nothing attenuating inside, the
  closed-form per-shape path runs unchanged (`has_attenuating_interior()` is cached).
- **Trace vs self-attenuating** are not separate shapes: efficiency depends only on the
  spatial distribution and the medium; InterSpec's per-volume/mass concentrations are
  downstream activity normalization.

## Environment Geometry (air + ground/walls) — offline scatter tier

A **third** geometry tier, `EnvironmentGeometry` (`src/geometry/EnvironmentGeometry.{h,cpp}`),
models the space *between* the source and detector tiers: homogeneous air plus N material
slabs. **0 slabs** = open air, **1 slab** = a field measurement over ground (source and
detector at any heights), **6 slabs** = a room (`::room(...)`). An infinite ground
half-space is the degenerate room, so both share all tracing and scoring code.

**This tier is an offline truth generator; it is NOT wired into `compute()`.** Routing a
full analog environment through the per-event efficiency loop would add cost *and* an
FEP-coupling hazard for a peak effect the analytic air attenuation already covers. Its
runtime product is the fitted environment-scatter payload (Full-spectrum approximation).

- **Next-event (point-detector) estimator** (`src/transport/NextEventEstimator.{h,cpp}`):
  scatter-in is wide-angle, so cone-biased brute MC cannot converge the scattered flux at a
  small far detector. The estimator scores, at **every** vertex of an analog
  multiple-scatter walk, the attenuated probability of reaching the detector. Multiple
  scatter is mandatory (air at 50 m / 662 keV has τ ≈ 0.46, double/single ≈ 23%; slab
  albedo is intrinsically multi-collision).
- **The product is a LINEAR incident-fluence tally** φ(E′, μ_in) (`FluenceTally`), not an
  efficiency. Peak transmission T(E0) and the scattered continuum are mutually exclusive
  outcomes of one primary and are never co-tallied (that would reintroduce the
  non-linearity hazard of Event Biasing). The tier delivers incident flux only, split by
  entry surface via `io/EnvScatterKernel.h` (face / side / back, from exact bin-averaged
  projected areas); detector response is applied downstream.
- **Bounded estimator**: the detector is a *sphere* and the score is expected **track
  length** through it (φ = L/V), with in-sphere attenuation `(1 − e^{−μc})/μ`; the
  constructor throws if a slab reaches within R of the detector. The attenuation leg
  follows the **sampled** direction, not the direction to the sphere centre (Previously
  tried lists the forms that biased it).
- **Uncertainty**: `FluenceTally::total_sigma` is exact from per-history totals. Summing
  per-bin variances in quadrature under-states it ×1.485 in production (one history scores
  many bins with positively correlated magnitudes); `surface_flux_sigma` still does that
  (TODO.md).
- **Internal validation**: `NextEventEstimator::run_brute` is an independent analog truth
  (real track length in the sphere) sharing only geometry and physics kernels:
  high-acceptance agreement 0.01%; the six-surface concrete room
  (`estimator_matches_brute_force_in_a_room`) agrees on the total (z ≤ 1.3) and on the
  multi-slab chain fraction (0.341 vs 0.346); the annihilation channel is gated (z −0.72 at
  2614 keV) and the brems channel at the total and channel level
  (`estimator_matches_brute_force_with_brems`). These live in `tests/test_env_scatter.cpp`
  behind `-DCEELO_RUN_MC_TESTS=ON` (~12 min). Shared-mode errors in the common kernels are
  invisible here by construction — that is what GEANT4 is for.
- **GEANT4 validation — against `G4EmLivermorePhysics`, the harness default, NOT
  option4** (see the harness's reference-physics note for how the two lists differ).
  `ceelo_g4val --sphere-flux R` scores the same observable (φ = L/V in an air sphere, all
  generations, direct/scattered split, creator-process origins);
  `tools/full_spectrum/env_scatter_export.cpp` writes the matching GDML and macro. Do NOT
  use `--fluence` here — it is primary-only and first-crystal-entry, so it drops the
  secondaries that are the signal. Normalization is confirmed by the analytic direct term
  (GEANT4/analytic = 1.005, 1.004, 0.996). Scattered flux, 60 m world, CeeLo at 4M
  histories:

  | Config | GEANT4 relative to CeeLo | z |
  |---|---|---|
  | field/soil 662 keV, 5 m, h = 1/1 m | +0.30% | +0.16 |
  | field/soil 662 keV, h = 0.3/1 m (72M G4 events) | +1.54% | +1.13 |
  | concrete room 662 keV, 2 m | −0.57% | −1.72 |
  | field/soil 122 keV | +2.37% | +1.30 |
  | field/soil 2614 keV | +0.24% | +0.11 |

  The brems channel alone at 2614 keV: 4.40e-9 ± 0.09e-9 (CeeLo) vs 4.49e-9 ± 0.27e-9
  (GEANT4), z −0.3. The room's agreement externally confirms the multi-surface fraction
  below. Only **integrated** scattered flux is compared; the energy/angle shape is not yet
  (TODO.md).
- **THE GDML WORLD SIZE IS PART OF THE COMPARISON.** GEANT4 kills tracks leaving the world
  while `EnvironmentGeometry::field()` is an infinite half-space, so a small world charges
  CeeLo for flux the reference cannot produce — a one-directional bias indistinguishable
  from a physics excess (a 15 m world read 3–5 points low at 122–662 keV; contributions
  fall off ~r⁻³ per unit radius, so the binding scale is geometric, not the air mean free
  path). **State the world size on both sides of any open-geometry comparison.** The
  exporter default is `max(12·distance, 6000 cm)`; `EnvScatterConfig::world_half_size_cm`
  (default 0 = unbounded, honoured by `run_brute` too) matches a small world deliberately.
- **Secondaries**: brems from every walk electron (Compton recoil, photoelectron, both pair
  electrons) via `ElectronCsda::sample_brems_emissions_in_material` (true along-walk
  emission positions — metres matter in air), scored next-event through the
  sampler-faithful 2BS density (`physics/BremsAngular.{h,cpp}`) about the electron's
  diffused direction, then walked with a `FluxOrigin` tag mirroring GEANT4's creator
  semantics. `include_brems` (default on), `brems_electron_floor_keV = 150` (the channel is
  4.26% / 0.29% / 0.04% of scattered flux at 2614 / 662 / 122 keV), `walk_brems` A/B knob;
  cost ≤2%. Both annihilation photons are walked from where the positron stops (positron
  ranges are metres in air; 40% of pair vertices at 2614 keV are in air); annihilation in
  flight is not modelled (~2–3% of positrons).
- **Discrete LINE channels** (`FluenceTally::add_line`, `LineOrigin`): flux at an exact
  energy stays out of the continuum grid, because a ~10%-wide log bin would make a
  photopeak indistinguishable from continuum downstream. **ElasticE0** (chains with only
  coherent scatters — 3.55% of scattered flux at 60 keV, half the content of its own bin;
  0.47% at 122 keV) and **Annihilation** (the 511 pair). The boundary is exact (zero
  inelastic events), and `total_mean()` still includes the lines.
- **Supporting pieces**: `RayleighAngular.{h,cpp}` evaluates the Rayleigh angular density,
  which needs F(x,Z) at an arbitrary angle, so the runtime photon data carry F from
  **EPDL2023 MF=27 MT=502** (`RUNTIME_ANGULAR_SPECS` in
  `tools/prepare_cross_sections/generate_epics_photon.py`, read through
  `CrossSectionData::form_factor_F()`); `AngularNorms.{h,cpp}` holds per-Z log-E norm
  tables (an exact-key cache is useless once energy varies at every vertex, which is also
  why `ComptonBiasConfig::max_biased_vertices` is 1); `GaussLegendre.h` the shared GL32 +
  s = sin(θ/2) substitution.
- **Drivers**: `tools/full_spectrum/env_kernel_bake.cpp` (`--mode air|field|room`,
  `--brute N` for the analog cross-check, `--fit-block` for the `CEEFIT-BLOCK-v1` blocks the
  fitted payload is built from) and `env_kernel_merge` (reassembles sharded bakes; seeds
  derive from the global node index, so sharding cannot change values).
- **What the environment does**: at 662 keV, 5 m, 1 m above soil, scattered flux is
  **0.675× direct per entering photon / 0.450× per scalar fluence** (always say which basis;
  the fluence one compares to published number albedo, 0.3–0.5 for soil/concrete), and
  **85% comes off the ground**. In a small concrete room it reaches 1.76× direct and
  **18–34% comes from chains touching two or more distinct surfaces**, growing as the room
  shrinks — so superposing single-surface kernels under-predicts by that much. Single
  scatter is 55.5% / 43.2% / 39.2% of scattered flux at 122 / 662 / 2614 keV and 17.2% in a
  room, which bounds any "deterministic single + tabulated multiple" basis.
- **It is a CONTINUUM effect, not a peak effect**: 0.3% of the scattered flux lands within
  ±10 keV of the 662 keV line and 97% is below 500 keV (albedo spectra peak at 150–300 keV).
  The peak is handled by the analytic transmission; what the environment adds is
  low-energy continuum (MDA, continuum under peaks, template confusion).
- **Why not GADRAS's scatter model**: against GEANT4 for air between a point source and a
  bare 3"×3" NaI (59–2614 keV, 5 cm–50 m), InterSpec's `GadrasShieldScatter`
  *transmission* agrees to <1% at ≥122 keV, but its *scatter continuum* is ~1.7–2× low
  near/mid-field with a ~20% sensitivity to the air extent it cannot parameterize; its
  environment term is one empirical clutter knob. Cascades need no special air handling for
  the peak (a sum peak scales as T(E₁)·T(E₂)); scattered-plus-full-energy coincidences are
  open (TODO.md).

## Full-spectrum approximation

**Purpose.** The Monte Carlo engine answers one configuration in seconds; InterSpec needs
pulse-height spectra and total efficiency for arbitrary (source, shield, distance,
environment) combinations in milliseconds, differentiable in the shield parameters, so
that activity and shielding can be fitted and cascade summing can use a total efficiency
that includes shield and environment scatter. The approximation composes separately
fitted operators, each a compact coefficient payload fitted to labelled CeeLo Monte Carlo
truth:

```
photon lines ──► SourceGeomFluxModel ──────────────┐   escaping photons per decay:
e-/e+ lines  ──► ChargedParticleSourceFluxModel ───┤   SpectralMeasure (atoms + binned density)
                                                   ▼
      [environment: EnvFitPayloadEvaluator -> scattered φ(E, μ) at the detector — not yet composed]
                                                   ▼
          geometric acceptance × finite-distance factor × DetectorPulseResponseModel(E_in, cos θ)
                                                   ▼
          Gaussian resolution, FWHM² = a0 + a1·E + a2·E²  ──►  channel counts
```

Code: `src/full_spectrum/SpectralMeasure.h` (the measure),
`src/full_spectrum/FullSpectrumModel.{h,cpp}` (the source, charged-particle and detector
operators and the composer `FastSpectrumModel`), `src/full_spectrum/ArrivalResponse.{h,cpp}`
(arrival-averaged response and the fold table), `src/io/EnvFitPayload.{h,cpp}`
(environment runtime), `src/io/EnvFitArchive.{h,cpp}` and `src/io/EnvFitBlockBuilder.{h,cpp}`
(environment transport blocks); bake and fitting tools in `tools/full_spectrum/`, the
environment campaign package in `tools/full_spectrum/env_scatter_fit/`. Examples:
`examples/full_spectrum_nuclides.cpp` (nuclide histograms),
`examples/full_spectrum_ceres_fit.cpp` (activity/shield fits).

### Representation rules

- **Atoms are atoms.** Every line-like outcome — uncollided and Rayleigh-only lines,
  511 keV, escape peaks, Ge K-escape and absorbed K lines, the full-energy peak, the
  environment's ElasticE0 and 511 lines — is a discrete `SpectralAtom`, never smeared into a
  continuum basis. Continua are integrals of positive densities over caller-requested bins
  (`SpectralBin`); no stored histogram is evaluated.
- **Rate × normalized shape** per component, positive by construction: log-yields, softmax
  or nonnegative archetype mixtures, nonnegative B-spline coefficients.
- **Steep dependence goes in a mechanistic base; the fitted field carries a smooth
  residual.** The fit stores log(truth / base) and the evaluator multiplies the base back
  in, so a C++/Python disagreement about the base is a multiplicative bias no fit metric
  can see. Every base, feature and basis therefore has ONE implementation in C++; the Python
  fitters read them through helper binaries (`ceelo_source_cross_section_features
  --emit-component-bases`, `ceelo_rate_surface_anchor`), guarded by
  `check_base_lockstep.py`, `check_feature_lockstep.py` and `check_basis_lockstep.py` (the
  bake's build step runs all three).
- **Differentiable.** Evaluators are `Scalar`-templated with no Ceres dependency;
  `ceres::Jet` overloads participate through ADL. `tests/test_full_spectrum_ceres.cpp`
  gates Jacobians against central differences.
- **Labelled truth, compared by family.** Components are tallied by interaction-history
  label in the CeeLo campaigns and the GEANT4 harness (`--component-labels`). The two codes
  count Compton vertices differently (at 800 keV the single/multiple-Compton labels differ
  by −/+15% while their sum agrees to 0.17%), so compare component *families*.

### Source escape (`SourceGeomFluxModel`)

Photons escaping a centred point source in a homogeneous sphere, per emitted photon. The
shield is described by areal density (g/cm²), effective Z — a Gaussian cardinal blend of
element cross sections over Z = 1–92 (`z_blend_width` 0.42), so Z is continuous and
differentiable — and hydrogen mass fraction f_H, modelled as exact Z = 1 mixed with the
effective-Z remainder. f_H is an input, not a free fit parameter (for light hydrogenous
shields areal density, Z and f_H collapse to one electron-areal-density parameter).
Validity 10 keV (trained from 59 keV) to 12 MeV, Z 1–92, 0–60 g/cm², optical depth trained
0–8, f_H 0–0.25.

- **Uncollided** atom at E0: exp(−τ_total) (Rayleigh included when `enable_rayleigh`).
- **Rayleigh-only ("elastic") line** atom at E0:
  `elastic_escape_scale · e^{−τ_abs} (1 − e^{−τ_rs}) e^{−τ_cs}` — the Compton damping is the
  physical 1.0 and the scale is fitted (`fit_elastic_escape.py`).
- **Single and double Compton**, exact: quadrature over the KN × S(x,Z) angular law, the
  bound-electron (Doppler) energy distribution with CeeLo's subshell profiles, and
  interaction depth; Compton+Rayleigh mixtures are separate components.
- **Three or more Compton scatters** — `MultipleScatterSurrogatePayload`: a Chebyshev-in-τ
  conditional yield (frozen beyond the fitted τ maximum) and RBF-interpolated shapes over
  (log E, Z/92, log τ).
- **Secondaries** (`SourceSecondarySurrogatePayload`), six terminal components (and the
  shield x-rays below):
  `pair_511_atom`, `pair_annihilation_continuum`, `pair_annihilation_in_flight`,
  `photoelectron_brems`, `compton_electron_brems`, `pair_electron_brems`. Each is a
  mechanistic base (process onset × escape, `detail::source_component_bases`) × a fitted
  conditional yield (local-linear RBF over the 7-D process coordinate: log E, Z/92, log τ,
  parent-process optical depth, the descendant's own absorption depth, thickness in
  radiation lengths, f_H) × a normalized cubic-B-spline shape. The escape term is the depth
  integral `(e^{−τ} − e^{−κτ}) / (τ(κ − 1))`; for broad-spectrum children (brems, in-flight
  annihilation) it is **averaged over the child's own emission spectrum** (GL8 in ln k over
  [10 keV, k_max], removal PE + Compton + pair; `spectrum_averaged_escape_factor`). Beyond
  the deepest populated τ the conditional yield decays as exp(−κ_tail (τ − τ_last)), with
  `tail_decay_kappa` = 0.8 set by the bake (the struct default 0 disables it). Descendants
  are fitted terminally, not composed recursively. The brems onset is a smooth turn-on in
  the parent electron's energy, 1 − e^{−(T − 200)/650} (the walk radiates only above
  200 keV); the photoelectron's T is the photon's less the binding of the shell it left,
  element by element of the blend as the walk samples it (K, retained L, else none), so a
  K-shell photoelectron from 250 keV in lead (162 keV) does not radiate.
- **Shield K and L x-rays**, three more secondary components — `fluorescence_k_alpha_atom`,
  `fluorescence_k_beta_atom`, `fluorescence_l_atom`: the shield's own x-rays that escape
  unscattered, yield-only atoms at the element's lines (lines below 10 keV are neither
  emitted nor tallied). The base sums, line by line, the vacancies made at each depth of the
  primary's path times the line's escape from that depth (the spherical-chord average with
  PE + Compton + pair removal at the line, so a Rayleigh-only x-ray stays in the atom). The
  vacancies come from the primary's photoabsorption (K and L1–L3 by their photoelectric
  shares) and its Compton collisions (shells by occupancy), a once-Compton-scattered
  primary absorbed on its way out, absorbed annihilation photons, impact ionization by the
  primary's photo-, Compton and pair electrons, their bremsstrahlung absorbed in the shield,
  and, for the L lines, the shield's own K x-rays absorbed where they are born. Depth is
  GL8 segments refined geometrically toward the surface (an opaque line escapes from a skin
  1e-4 of the radius thick) and the escape angle is split at grazing. Each element of the
  effective-Z blend carries its own lines, its weight tapered C¹ to zero below 0.5% and
  renormalized, so an integer Z emits exactly its own lines. A family has one fitted
  conditional yield on the process coordinate (parent: the primary's photoabsorption;
  descendant: the family's intensity-weighted removal), spread over its lines by their
  bases. On a 198-cell span bank (Z 32–92, 59 keV–12 MeV, τ 0.2–8) truth/base reads 1.00–1.04
  at every τ just above the K edge (100 keV in W and Pb) and grows with depth from ~660 keV
  up (Pb at 1.5 MeV: 1.3 at τ 0.2 to 3.3 at τ 2): photons Compton-scattered more than once
  and then absorbed make K vacancies the base does not count. Median span 2.0 / 1.8 / 1.25
  (Kα / Kβ / L), slope +0.29 / +0.31 / +0.08 per unit τ; the K families carry their own
  slope bounds in `check_kernel_span.py`. An element's lines are compressed once per family: lines under 0.5% of the
  family's yield fold into the nearest kept line, and lines within 0.05 keV merge.
  `predict` merges the atoms across source lines, which share the shield's line energies.
  Compton-scattered x-rays are tallied by the campaign (`fluorescence_continuum`) and
  reported against the atoms, but not modelled (TODO.md).
- **Photonuclear** (`PhotonuclearSourcePayload`): mechanistic (1 − e^{−τ_pn}) with a
  closed-form escape factor × a fitted (E, Z, τ) conditional yield and shape over a
  natural-element GEANT4 photonuclear table. **Off by default** (`enable_photonuclear =
  false`) and baked only on request (`PHOTONUCLEAR=1`); it matters only above ~8 MeV.

### Charged-particle source (`ChargedParticleSourceFluxModel`)

Photons escaping a centred monoenergetic electron or positron in an elemental-equivalent
sphere (kinetic energy 10 keV–12 MeV, Z 1–92, 0–1000 g/cm²; hydrogen via Bragg-additive
ranges and linear X0⁻¹). Thirteen components: `electron_brems`, `positron_brems`,
`positron_containment`, `positron_511_atom`, `positron_annihilation_continuum`,
`positron_annihilation_in_flight_direct`, `positron_annihilation_in_flight_continuum`, and
the shield's unscattered K-alpha, K-beta and L x-rays for each particle
(`{electron,positron}_fluorescence_{k_alpha,k_beta,l}`). Each is a mechanistic yield base ×
exp(a tensor clamped cubic B-spline log-yield field over log-KE warp, Z/92 and warped range
ratio) × a rank-R (≤ 16) softmax mixture of nonnegative archetype shapes (cubic B-splines in
log photon energy, under an absorption envelope); the x-rays are yield-only atoms at the
element's lines. Coordinates freeze-clamp to the trained domain, so outside it the
mechanistic base carries the trend. `predict` feeds continuous beta branches as a fixed
positive quadrature of `ChargedParticleLine`s and conversion electrons as lines. Fitted by
`fit_charged_particle_source_structured.py` (yield fields) and
`fit_charged_particle_joint_shape.py` (one alternating minimization of the deployed shape
error; settings in `tools/full_spectrum/bake/config/`).

The x-ray base is a + contained × b (`detail::charged_fluorescence_bases`): a from impact
ionization along the particle's path and from its bremsstrahlung absorbed in the shield, b
from the two annihilation photons, per unit positron containment. The path runs radially
out to the stop radius (the smooth min(1, R/t) of the annihilation atom), and the
impact-ionization vacancies sit on it. A photon born on the path is followed along its own
ray (8 sites × 16 directions), and its line escapes from where it is absorbed: at depth d
below the plane tangent at the exit, as (E₂(μ_x d) + E₂(μ_x (2t − d)))/2 — exact in a thin
shield, the half-space in a thick one, and closed-form along each ray. In a shield thicker
than the range the x-rays that escape come from bremsstrahlung absorbed near the surface,
far from where it was born. The bremsstrahlung count shares the brems bases' Kramers / X0
normalization and straight path, both low for slow electrons in high-Z shields, so on the
smoke bank truth/base is 1–6 except for Pb L lines from a 100 keV electron (~45).

### Detector response (`DetectorPulseResponseModel`)

Deposited-energy measure per photon incident on the detector **envelope** (endcap
included), for a parallel beam at (E_in, cos θ). Atoms: `zero_deposit`, `full_energy`,
`pair_full_energy`, `single_escape`, `double_escape`, `annihilation_511`, Ge
`k_alpha/beta_escape` (at E_in − 9.833 / 10.943 keV, after any number of Compton scatters
and from a photoelectric or, with spectrum fidelity, a Compton K vacancy) and
`k_alpha/beta_absorbed` (fixed 9.833 / 10.943 keV deposits: the vacancy is made in a dead
layer and only its K x-ray reaches the active volume). The banks run with spectrum fidelity
on. Continua: `photoelectric_escape` (one linear basis) and the
Compton and pair-history continua, each expanded in **independent cubic B-spline blocks
between moving physical landmarks** — [0, Compton edge] and [edge, E_in] for single and
multiple Compton; sorted {0, 511, E_in − 1022, E_in − 511, E_in} for the pair continua —
with a fixed coefficient count per block so a coefficient keeps its meaning as E_in moves.
Pair outcomes are hierarchical: `pair_fraction` holds the odds of pair history among
detected pulses and pair outcomes are normalized within that population. Component yields
and shapes are normalized-RBF fields over (log E, θ/π, passive-layer landmark
1/(1 + (E/E_L)³)), with a threshold-expanded energy coordinate for pair-history components.
A yield field's stored node values are coefficients solved so the normalized kernel
average reproduces the training values at the nodes: 16 Landweber steps v ← v + (y − S v),
S the row-normalized kernel matrix at the final bandwidth (`interpolate_yield_nodes`);
S's eigenvalues lie in [0, 1], so the field at the nodes is the data filtered by
1 − (1 − λ)^17 ≤ 1 (no noise amplification), and the evaluator is unchanged. Shapes are
plain kernel averages. The detected rate is a global Bernstein log-odds surface (32 energy
× 24 angle terms, angle warp x(x² + ε²)^((q−1)/2), q = 0.5, ε = 0.03) plus a deterministic
surface-entry transmission anchor below ~150 keV (`rate_surface_anchor_log_odds`; no fitted
parameters, built from the geometry and the mass-attenuation tables). K-line atoms are
gated off above per-component support limits, because a normalized-RBF yield levels off
instead of decaying outside its node cloud.

**A payload describes one detector, and carries it.** The `CEELO_DETECTOR_RESPONSE`
payload (and the `detector_*` records of `CEELO_FULLSPEC`) embeds the detector's label and
its `GeometryDescriptor` text; the reader rebuilds from it the crystal, the anchor's
geometry and the envelope, so nothing geometric is stored twice. The anchor
(`set_rate_surface_anchor_geometry`) follows `trace_ray()`'s shell stacking: front, side
and back passive stacks (a layer covers the back only where its `z_end` reaches past the
shell inside it; a collimator joins the side stack), the dead layers, and the bore's
opening in the back face, which backward photons enter onto active crystal. P0 is per
photon on the **envelope** (`outer_bounding_radius`, `outer_z_extent`), the campaign's
normalization, so its face, bore and side terms are divided by the envelope's projected
area πR_e²|cos θ| + 2R_eL_e sin θ; `predict` uses the same area for its acceptance.

The landmark E_L is the energy at which the front entry (stack plus front dead layer,
normal incidence) has optical depth `kPassiveLandmarkOpticalDepth` = 0.770421,
`hpge_coax`'s depth at the 55 keV its feature map was selected with
(`passive_landmark_keV`; 62.2 keV for the Micro-Detective housing). It is stored at full
precision, because the fitted node features were computed with it. The fitter reads the
detector from the bank manifests (`run_detector_response_campaign.py` records
`--print-geometry`) and asks `ceelo_rate_surface_anchor --geometry FILE --info` for the
landmark, K-line element and K edge, so it carries no detector constants of its own.

The detector bake (`DETECTOR=N tools/full_spectrum/bake/detector.sh`) fits the built-in
`hpge_coax` (8.0 cm diameter × 8.0 cm, 0.5 cm-radius 5.5 cm-deep bore, dead layers 0.05 /
0.08 / 0.08 cm front / side / back, 1 mm Al endcap) or any `bake/config/<N>.geometry`;
`detective_x` is a placeholder (a Detective-X-like crystal in the Micro-Detective
housing). `nai3x3_payload()` (the default) and `coaxial_hpge_payload()` are hand-set
mechanistic pilots, not fits. Truth comes from `ceelo_detector_response_campaign` with
`--detector-brems --detector-fluorescence --detector-spectrum-fidelity` (the last is
`enable_spectrum_fidelity`: Compton K-vacancy x-rays and the Detailed electron walk);
the campaign's own defaults leave all three off, which is wrong for a payload (Previously
tried). The bake places its banks per detector
(`detector_bank_energies.py`): energies below 100 keV are written for `hpge_coax` and
scaled by E_L / 55 keV, so the `lowE` turn-on bank (0.64–1.32 E_L) and the low audit and
arrival energies sit at the same place in any detector's turn-on; every held-out energy is
kept ≥ 2% from every training energy. The placed recipe is recorded, and a resume under a
different recipe is refused.

**Arrival response (`ArrivalResponse.h`).** Scattered photons arrive from every side, so a
scattered field is folded with the fixed-angle response averaged over arrival direction,
per envelope photon weighted by the projected area A(μ) = πR²|μ| + 2RL√(1−μ²). An arrival
band is an isotropic field restricted to direction cosines [μ_lo, μ_hi]: iso [−1, 1],
front [0, 1], back [−1, 0], side [−0.3, 0.3]. Both sides of the average live in the library.
`ArrivalBeam` is the Monte Carlo side (the campaign's `--incidence`): μ and azimuth uniform
in the band, the origin uniform on the bounding sphere's disk, both redrawn on a miss, so
accepted photons are distributed as A(μ) and a tally per accepted photon is per envelope
photon. `arrival_quadrature` / `evaluate_arrival` are the model side: Gauss–Legendre in θ
with 16 nodes per hemisphere (A(μ) has a kink at μ = 0 that a spanning rule resolves
slowly; per hemisphere GL-16 is within 1e-5 of converged). The bake's held-out `arrival`
bank (15 energies, 0.75–1.64 E_L and 110–850 keV, × the four classes, 4M events per cell)
is scored by `score_arrival_bank.py` and gated per class and energy band: detected yield,
continuum TV against its MC noise floor, the valley between the Compton edge and the line,
and the envelope area. Arrival cells carry `cell_kind detector_arrival`, which the fitter
refuses.

**A point source at a finite distance (`FiniteDistanceDepth`).** The response is per
envelope photon of a parallel beam, so the far-field acceptance (projected area over 4πd²)
overstates a point source at d: its rays diverge, and photons entering near the rim leave
through the side first. Each payload carries, per energy, a virtual disk — depth δ and
radius R — for the atoms that track the incident energy and another for everything
detected, and `predict` multiplies the far-field response by 2d²/(s(s + x)), x = d + δ,
s = √(x² + R²). The bake's `finite_distance` step (`ceelo_detector_finite_distance`)
measures them by point-source MC of the payload's own geometry on the axis (21 energies ×
13 distances, 5 cm–10 m, 1M events per point, up to 16× where noisy), with distances from
the geometry's reference point as everywhere in CeeLo (`DetectorResponse::query_position`).
It fits (δ, R) by grid search so that ε(d)·s(s + x) is constant on the ladder, keeping a
fit only if every point has ≤2% noise and the spread is ≤2% (else the nearest kept energy
stands in, below ~30 keV on axis). For `hpge_coax` the full-energy disk goes from δ 0.06 cm,
R 4.7 cm at 50 keV to a point at δ 3.4 cm above 3 MeV, with fit spread 0.04–0.25%. The
depths are applied for incidence cos θ ≥ 0.5 (`kFiniteDistanceMinCos`) and the far field
below.

**Folding an incident field (`ArrivalResponseTable`).** The same payload, tabulated for
folding a scattered field φ(E, μ) by one contraction: per arrival-cosine node (the same 16
GL nodes per hemisphere) and incident-energy node (log-spaced, 2%), the continuum and the
fixed-energy atoms (511 keV, the absorbed K lines) on the caller's deposit grid, and the
yields of the atoms that track the incident energy. An `IncidentField<Scalar>` holds
fluence per (energy cell, μ bin) and per line; only those values carry `Scalar`, so a fit's
Jets reach the fold through the field alone. Energy mass is spread over the nodes by a hat
basis; in μ the envelope's exact projected area is integrated against the response
interpolated by the polynomial in θ through each hemisphere's nodes (barycentric
Lagrange); tracking atoms land on their own interval [E_lo − offset, E_hi − offset], so a
smooth continuum never renders as a comb of node lines. Against a converged direct fold of
air-and-ground kernels on `hpge_coax` (4096 channels to 3 MeV): integral within 0.012%
(Cs-137, Co-57, Th-228 at 2614.5 keV, Am-241), TV 0.06–0.08% raw, worst bin within 1.1% (a
midpoint-rule fold per kernel cell is 2.0–2.7% TV). Release, M1 Pro: tables build in
0.15–0.66 s to 122–662 keV and 4.7 s to 2.6 MeV (0.5–24 MB); a fold of a 384 × 32-cell field
takes 4–9 ms, 6–18 ms with 8-derivative Jets.

### Environment scatter (`EnvFitPayloadEvaluator`)

A fit-first campaign (`tools/full_spectrum/env_scatter_fit/`, whose README is the
contract): a deterministic Sobol design of nodes over (E0, distance, source height,
detector height, extent) with disjoint training / validation / final roles (256 / 64 / 64
in production); each node is transported by the next-event estimator in `CEEFIT-BLOCK-v1`
blocks of 262,144 histories (32 replicas) and ingested into a `CEEFIT-HDF5-v1` archive that
keeps per-history sufficient statistics. Each continuum component is fitted to ŝ = raw
fluence / analytic direct fluence as a positive log-yield and softmax geometry weights on a
4-D clamped tensor cubic basis, times joint energy–angle archetypes (nonnegative mixtures of
unit-integral cubic B-spline PDFs on a **registered** energy axis — log blocks split at the
backscatter-edge landmark, so the edge sits at a fixed coordinate for every E0 — and
t = (μ + 1)/2); ElasticE0 and Annihilation-511 are exact lines with angular models. Rank
and level are selected on validation nodes; adaptive rounds allocate more histories (a
continued node at most doubles per round) or nodes by predicted validation loss per CPU,
and the campaign stops after two consecutive passing rounds. Blocks 1–4 of each node stay
in the archive as raw records (the uncertainty calibration pairs them); later blocks fold
on ingest into the node's running sums, kept as the fold's exact state (every Kahan total
and compensation, the integer sums, a checksummed ledger row per block), which reproduces
the fold over all blocks bit for bit (the field/soil campaign: about 25 GB, at most ~70 GB,
instead of ~456 GB). Each round's seal re-verifies only the nodes it changed (a
`.verified-nodes.json` ledger); publication requires a full seal. The runtime
`CEEFIT-PAYLOAD-v1` is dependency-free and checksummed; `EnvFitPayloadEvaluator` takes the
caller's analytic direct fluence T(E0)/(4πd²), returns ŝ and incident-fluence maps per
component and line, and rejects extrapolation outside the design box and
`|h_src − h_det| ≥ d`. A Rayleigh-on environment payload pairs with the
`AnalyticCollimated` direct term. Production domain: E0 20–3000 keV, d 0.25–25 m, h_src
0–4 m, h_det 0.3–4 m, 384 energy × 64 angle bins, Doppler off; tiers field/soil (default),
field/concrete and open air.

### Payload formats

Every full-spectrum format and dataset is **version 1**, with no backward-compatible
readers, until the first release: a format change edits the one reader and the payloads
are regenerated. There is no provenance machinery (no recorded code revision, dirty-tree or
physics hashes); a bake runs from one code version and the commit that adds its payload
says which machine and dates. Integrity checks that feed identities or seeds (CEEFIT
checksums, node IDs, model-content hashes) stay.

| Header | Content | Loaded by |
|---|---|---|
| `CEELO_SOURCE_ESCAPE 1` | source-escape model: validity, scalars, elastic line, multiple-scatter surrogate, secondary surrogate (with the shield x-ray families; a component without nodes emits nothing), photonuclear section (each section `model_version 1`) | `SourceGeomFluxModel::deserialize` |
| `CEELO_CHARGED_PARTICLE_SOURCE 1` | thirteen spline components, each required (layout in `FullSpectrumModel.h`) | `ChargedParticleSourceFluxModel::deserialize` |
| `CEELO_DETECTOR_RESPONSE 1` | detector geometry, RBF components, block sizes, rate surface and anchor, K supports, finite-distance depths | `DetectorPulseResponseModel::deserialize` |
| `CEELO_FULLSPEC 1` | a composed `FullSpectrumPayload` (source section + detector + distance, incidence, quadrature) | `FastSpectrumModel::deserialize` |
| `CEEFIT-PAYLOAD-v1` (binary) | environment runtime model and its model inventory | `load_env_fit_payload`, `EnvFitPayloadEvaluator` (`src/io/EnvFitPayload.h`) |
| `CEEFIT-BLOCK-v1` (binary) | one transport block (raw fluence moments, replicas, stopping vector) | `src/io/EnvFitArchive.h`; the Python archive |
| `CEEFIT-HDF5-v1` | the environment scientific archive: node aggregates, raw records of each node's first four blocks, per-node fold states | `env_scatter_fit` only |
| `CEEK` binary 1 / text `# format: 1` | tabulated environment kernel grid (GEANT4 comparisons, cross-checks) | `load_env_kernel` (`src/io/EnvScatterKernel.h`) |

### Regeneration and acceptance

Payloads are regenerated from this source tree alone by the drivers in
`tools/full_spectrum/bake/`, one per family, following
[tools/full_spectrum/README.md](tools/full_spectrum/README.md). A driver generates the
Monte Carlo truth banks (one gzip-compressed CSV per CeeLo cell), fits, scores the
candidate on held-out banks and independent checks, applies the family's gates
(`check_detector_gates.py`, `check_source_escape_gates.py`,
`check_charged_particle_gates.py`; the environment campaign's calibration, acceptance and
publication reports) and copies an accepted payload to `data/full_spectrum/<family>/`. Gate
bounds are documented bars or the last accepted payload's value with headroom,
max(1.25 × last, last + 0.5 percentage points), so a statistically equivalent re-bake
passes and a regression does not. `SMOKE=1` runs every step on tiny banks.
`data/full_spectrum/` is git-ignored; only the final payload of a bake whose gates all
passed is committed (`git add -f`), never banks, intermediate fits or smoke output. Bakes
run natively (macOS with Homebrew is the documented platform; Linux works). A default
source-escape bake needs no GEANT4: the photonuclear term (`PHOTONUCLEAR=1`, ~450 CPU-hours
of GEANT4 at 1.1–4.0 ms per event) and the GEANT4 leg of the 150-cell acceptance grid
(`GEANT4_CHECK=1`, a few CPU-hours) are opt-in.

### Last measured accuracy

Headline numbers from the most recent bake of each family (failure detail is in TODO.md).
Where a gate has no documented bar, its bound is the value below with headroom. The
source-escape and charged-particle numbers are from the Nov 3 2026 re-bake (shield x-ray
components; fluorescence, impact ionization, L lines and EADL K energies in the truth); the
detector paragraph is the last fidelity-off bake, with the first fidelity-on bake after it.

**Detector (`hpge_coax`), validation bank (88 cells, 84 scored):** detected yield error
median 0.10%, p95 0.40%; pulse-shape total variation median 0.95%, p95 2.36%; full-energy
atom median 0.28%, **p95 2.29% — failing the 1.5% bar** (0.23% / 1.11% below 3 MeV).
4M-event spot checks at normal incidence: +0.02% (122 keV), −0.16% (661.657 keV, z ≈ −2),
+0.24% (1332.492 keV, z ≈ +2.4). Audit banks: full-energy median 0.46% / 0.35%, p95 2.89% /
2.10%. Arrival bank (15 energies × 4 classes): detected yield median ≤ 0.16% in every class
and band (front 40–100 keV, 0.163%, misses its 0.15% bar), continuum TV 1.2–2.6× its MC
noise floor, valley density median 1.3–2.6% (worst cell 9.0%), envelope area |z| ≤ 2.8.
Every detector baked so far fails the full-energy p95 gate (`detective_x_fit` 2.28%).
~2.2 ms per 4096-bin response; ~2.5 MB. The CeeLo truth itself vs GEANT4 (matched geometry
and physics, 2M events, per envelope photon): lines −0.16% (z −0.7) at 800 keV, −0.25%
(z −0.9) at 2614.5 keV, +1.44% (z +4.3) at 5000 keV; detected total within 0.2% (measured
before the crystal-walk straggling changes).

The first detector bakes with spectrum fidelity in their truth (Nov 3: Detailed walk,
Compton K vacancies, impact ionization) fail gates set on fidelity-off truth, in the same
two families: `hpge_coax` full-energy p95 2.93% (median 0.33%), arrival 40–100 keV back /
side detected-yield median 0.169% / 0.177%, back valley-density max 11.2%, side
detected-yield excess max 0.86%; `detective_x` full-energy p95 3.53%, detected-yield p95
0.53%, five arrival medians at 40–400 keV (0.15–0.24% yield, 3.2–4.9% valley). Spot checks
pass (`hpge_coax` within 0.04 / 0.37 / 0.10% at 122 / 661.657 / 1332.492 keV). Neither payload is
published.

**Source escape (gates pass):** 150-cell acceptance grid (59 keV–3 MeV, collided photons
only, 104 cells above the leakage floor), |median| 0.74%, p95 8.28% from CeeLo truth (the
opt-in GEANT4 leg last read 0.98% / 8.43%). 75-cell off-grid-τ holdout (τ 0.35–6, 4M events
per cell): leakage median 0.35%, p95 4.42%; deep-τ bremsstrahlung median 9.2%, p95 67%;
shield x-ray atoms median 3.05%, p95 35% (43 cells; their Compton-scattered continuum, not
modelled, is 1.2% of them); spurious yield past the trained τ ≤ 1.1e-7 photons per primary.
Secondary per-component yield error: 2.68% on-grid, 3.54% off-grid τ, 5.04% held-out element
(Z = 13, 35, 64). Multiple-scatter surrogate test bank: yield error median 1.16%, p95 8.75%;
shape TV median 2.34%, p95 8.40%. Photonuclear, when baked: 0.59% median at trained elements, 13.4% off-grid, 27.2%
held-out element (≤1.9% of escaping photons). CeeLo truth vs GEANT4: ≤0.2% on 25 cells of
{C, Al, Fe, Cd, Pb} × {0, 1, 5, 15, 50} g/cm² at 3 MeV, and ≤0.1% (z ≤ 0.3) for
hydrogenous shields (water at 300 / 3000 keV, polyethylene at 662 keV).

**Charged-particle source (gates pass), whole-Z holdout planes (Z = 35, 64), median / p95:**
shape TV electron_brems 2.68 / 14.4%, positron_brems 1.32 / 6.3%, annihilation_continuum
1.64 / 5.4%, in_flight_direct 1.05 / 5.8%, in_flight_continuum 3.6 / 8.7%, positron total
0.14 / 1.3%; yield error electron_brems 1.43 / 8.8%, positron_brems 1.68 / 12.2%; shield
K x-rays (yield error) 4.7–7.4% median on the holdout planes and 2.8–5.3% on the off-grid
final test, p95 13–67% (low-count cells). Off-grid electron banks, total-spectrum TV:
1.86 / 7.3% (validation), 1.39 / 6.4% (final test, Z = 16, 47, 79). C++ vs Python parity
1.2e-3. ~0.8 MB; 2.7–3.4 ms per 1024-bin call.

**Environment scatter:** the field/soil baseline has run (384 nodes × 4 blocks, 13.3 GB);
the adaptive rounds have not, and publication is blocked (TODO.md). The transport truth
agrees with GEANT4 Livermore within |z| ≤ 1.72 on all five configurations (table above).

**Composition:** `predict` against 4M-event `EfficiencyCalculator` runs of `hpge_coax`, bare
point source on the axis, 122–2614 keV and 5 cm–10 m: full energy −0.19% to +0.57%, detected
−0.03% to +0.61% (+0.3–0.5% at 30 cm–3 m, z 10–25: the depth fit's grid and its 1M-event
points). The far-field acceptance alone would be +1.9% (122 keV) to +6.6% (2614 keV) off at
1 m and +69–168% at 5 cm. At 59.5 keV the full-energy peak reads +0.9 to +1.2% at every
distance (the campaign-vs-`EfficiencyCalculator` offset, TODO.md). Source at 60°
(cos θ = 0.5), 662 keV: −0.72% / −0.45% / +0.03% full energy at 30 cm / 1 m / 3 m.
`predict` costs 3.8–4.9 s CPU for 100 lines × 4096 channels, almost all in the source model.

### Limitations

Open items are in TODO.md's full-spectrum areas. Out of scope for now: non-convex rooms,
non-isotropic source leakage into the environment, and position/azimuth structure of the
incident field (the environment interface keeps energy and one arrival cosine).

### Previously tried

Measured dead ends, so they are not re-explored. Each line: what was tried, its measured
weakness, what replaced it.

*Source escape*
- **Single proxy-energy escape for bremsstrahlung children** (evaluated at E/(1 + 8E/mₑc²),
  removal PE + pair, softened by a κ-scale knob). κ lands at 30–100, where the closed form
  collapses to e^{−τ}/(τ(κ − 1)) and its τ-slope stops depending on κ (the factor spans
  97,600× from τ 0.2 to 8, identically for Pb at 3 MeV, Fe at 3 MeV and Pb at 662 keV).
  Replaced by the spectrum-averaged escape with PE + Compton + pair removal; the two changes
  are complementary — median conditional-yield span (photoelectron / Compton-electron /
  pair-electron brems) 76.7 / 113.7 / 340.8× with no escape, 7.9 / 16.6 / 26.0× with the
  proxy, **1.4 / 4.5 / 2.6×** with both; deep-τ brems holdout 21.1% → 6.8% median. GL16 =
  GL8 to three decimals; a Kramers (T − k)/k weight is worse than 1/k; a 20 keV floor is a
  near-constant the fitted yield absorbs.
- **Creation-only bases** (`onset × (1 − e^{−τ_parent})`, no escape): the fitted fields had
  to carry a 355–1528× span, leaving the brems channels −75% to −81% at τ = 6.
- **A bare exp(−τ₅₁₁/2) escape term for `pair_annihilation_continuum`**: span 248×,
  τ-slope −0.882, against 9.2× and +0.291 with the depth integral.
- **Kernel-mean (Nadaraya–Watson) interpolation on the τ grid {0.2, 1, 4, 8}**: off-grid-τ
  yield error 43.6%. Densifying to {0.2, 0.5, 1, 2, 4, 8} alone gives 40.4%, local-linear
  alone nothing on the sparse grid; both together 4.3%. Kernel means also bias trained cells
  when the bandwidth is comparable to the node spacing.
- **Flat extrapolation past the deepest τ node, with log sentinels for starved cells**:
  spurious escaping yield up to 3.2 photons per primary at τ = 50. Poisson-bound anchors
  plus the κ = 0.8 tail decay bring the worst to 6.6e-6. The shape absorption envelope
  cannot do this job: it preserves yield by construction.
- **Hand-set elastic-line constants** (scale 0.74, Compton damping 0.35): +25.2% median,
  +3367% at τ > 6. The free fit lands on damping 0.9994, i.e. the physical 1.0.
- **Effective Z alone for hydrogenous shields**: total escaping yield 8.2% median / 278% p95
  on a 672-cell pilot bank, against 1.25% / 147% with the f_H descriptor. Fitting f_H freely
  is degenerate for light shields (activity recovers, but T, Z and f_H collapse to one
  parameter), hence f_H as a fixed input.
- **Recursive composition of descendants**: 31.8% median |yield error| against 4.8% for
  terminal fits on 36 held-out pair-brems cells.
- **Bandwidth selection on held-out elements plus the off-grid bank**: both arms were
  dominated by on-node cells, and it spent the only between-node holdout on model
  selection. Replaced by a leave-one-interior-τ-node-out CV, which over-widens instead
  (TODO.md); on the 150-cell grid, the one arbiter neither selected on, the two tie.
- **Photonuclear across effective Z**: leave-one-element-out is 23.1% median (giant-dipole
  secondary multiplicity follows nuclear structure, not Z). Without the mechanistic escape
  factor, off-grid τ cells were +66% / +74% (with it +11.5% / +16%).
- **Duplicate C++/Python base arithmetic**: a knob that reached one implementation and not
  the other gave a 1.07× (τ 0.2) to 1.73× (3 MeV, Z 82, τ 6) bias no fit metric showed —
  hence the single implementation and the lockstep checks.
- **Explaining the brems shape floor in thick high-Z shields** (TODO.md): it is not
  interpolation (τ = 4 at a trained Z is a trained node and the worst cell), not basis
  resolution (the inverse-energy warp puts a knot at 88 keV and 13 of 30 below it), and not
  the absorption envelope (disabled; and `loss_factor = min(max_ratio, E0/E)^power − 1`
  saturates at `max_ratio` = 10, so at 3 MeV it is constant below 300 keV and has no K edge).
  A mechanistic transmission envelope from the escape kernel's closed form is worse than the
  fitted shape (20–300 keV TV 19.4% vs 5.2% at τ = 3, 20.2% vs 10.6% at τ = 4).

*Charged-particle source*
- **Scattered normalized-Gaussian RBF fields** (23,280 nodes, 12.9 MB): off-grid total
  yield error 5.46% median / 11.38% p95, against 1.05% / 3.58% for the tensor B-spline
  fields (21.5× smaller). Whole-Z holdouts sharing the training grid hid it.
- **A four-stage shape fit** (per-cell NNLS → NMF archetypes → NNLS weights → penalized
  log-odds regression): errors compounded and a 1e-3 smoothness weight was ~30× the data
  term (electron_brems TV 5.33% median). The joint alternating minimization of the deployed
  shape error gives 2.58%, and rank 16 wins on both brems components.

*Detector response*
- **Truth without detector bremsstrahlung and fluorescence** (the campaign defaults): the
  line yield is too high by 0.2–1.5% up to 1 MeV, 3.8% at 2 MeV, 17.8% at 5 MeV and 46.0% at
  10 MeV while the detected total moves ≤0.2%, so every yield gate passes. Against GEANT4
  the 5 MeV line read +30.1% before and +0.5% after matching the physics.
- **One nonnegative B-spline over a warped [0, E_in] coordinate**: cannot represent the
  Compton-edge step (deficit below the edge, an edge spike, a spurious maximum);
  curvature regularization makes it worse. Landmark blocks: edge-region χ²/ndf 1098 → 0.73
  with the same 48 coefficients, 38% faster. Blocks narrower than a histogram bin give
  collinear columns; the fitter zeroes blocks with fewer than four interior bin edges.
- **An unwarped 20 × 12 rate surface**: angular error (+5% at |cos θ| = 0.97); the angle
  warp and 32 × 24 terms took the audit median from 1.42% to 0.15–0.25%. The plain
  sign(x)|x|^q warp has an infinite derivative at cos θ = 0, hence the ε.
- **Fixing the sub-landmark rate by coordinates or capacity** (a second energy logistic, a
  μ-based landmark, 40–80 energy terms, 2–16× low-energy sample weights): nothing beyond
  selection noise. The transmission anchor took the 42.5–45 keV p95 from 7.16% to 1.13%.
- **An anchor without the bore**: it saw only the back dead layer at backward incidence,
  where at 30 keV nearly all detection is through the bore (P0 / truth 0.03–0.2, a
  log-odds correction the smooth angle basis could not follow).
- **The energy warp as the high-energy knob**: every value large enough to help gives back
  the low-energy gain. Ridge damping of the coefficients above 2 MeV (κ 100, cut 0.75) takes
  11 MeV from 0.396% to 0.215% median; the rest needs training energies near 12 MeV (now in
  `dense_angles`).
- **More `photoelectric_escape` basis terms** (16 → 48): ≤0.05 TV gained, overfitting above
  300 keV, up to 14% slower. The real defect was the absorbed dead-layer K lines (66.5% of
  that continuum at 42.5 keV), now atoms: 42.5 keV whole-spectrum TV 2.89% → 0.74% median.
- **One shared K-line support limit**: the four components lose support at 320 / 175 / 280 /
  150 keV, so one limit lets the shorter-lived ones extrapolate flat over most of a decade.
- **Earlier layouts** — features (log E, cos θ, sin θ) and a single aggregate pair continuum
  — were replaced by the θ/π + passive-layer features, the threshold-expanded pair
  coordinate and the hierarchical pair model.
- **Training values as the yield nodes** (a plain kernel average): biased wherever a field
  curves within a bandwidth or its nodes crowd on one side — the full-energy atom read up
  to +3% (z +32) at its own training cells and 1.58% / 16.6% median / p95 on validation;
  the node solve brought that to ~0.3% / ~2.3%, the spot checks from +1.6 / +2.7% to within
  0.25%, and the arrival classes at 100–1000 keV from 1.0–1.4% to 0.33–0.37% median.
- **Inverse-density node weights** (against the 16-angle `dense_angles` energies
  outweighing sparse neighbours): no gain with node correction (validation full-energy p95
  below 3 MeV 0.93% vs 1.03%).
- **Piecewise-linear interpolation in μ for the arrival table**: the integral of a 59.5 keV
  air-and-ground field came out +0.47% with 16 nodes per hemisphere (+0.13% with 32, +0.04%
  with 64). The polynomial in θ through the same 16 GL nodes gives +0.012%.
- **A one-parameter virtual point detector for `predict`** (δ free, the disk radius fixed
  at the crystal's, as `src/io/VirtualDepthFit.h` fits): −1.2% to +2.4% off at 10–30 cm; a
  point detector (d/(d + δ))² misses by up to +15% at 5 cm at 122 keV. With R free as well
  the spread is ≤ 0.25% over 5 cm–10 m.
- **The on-axis depths side-on** (cos θ = 0): −16.5% at 30 cm, against −1.7% with the
  far-field acceptance; they are applied for cos θ ≥ 0.5 only.
- **Selecting hyperparameters on an audit bank** overfits it (audit max 11.96%). Select on
  the validation bank with a rule fixed in advance; audit banks are report-only.
- **Point deposits for source-continuum quadrature bins** in `predict`: they render the
  smooth escaping continuum as nq artificial peaks that move with nq. Box deposits of the
  bin width converge (nq TV 1.81% → 0.26%).

*Environment scatter*
- **Brute-force MC with cone biasing** for scatter-in: does not converge at a small far
  detector; replaced by the next-event estimator. A larger scoring sphere does not help the
  estimator (at 50 m, R = 15 / 30 / 50 cm give σ 1.20 / 1.24 / 1.45%); it helps only
  `run_brute`.
- **Crossing counts with the mean-chord identity** (4V/S holds only for a
  parallel/isotropic field; near vertices clip the sphere at grazing incidence: +4.3%);
  **no in-sphere attenuation** (~0.05% at R = 5 cm, ~+0.5% at R = 50 cm); **toward-centre
  attenuation legs** (for a vertex buried in the ground the most vertical, shortest escape
  path — under-attenuates; invisible without a dense slab); **a per-bin-quadrature total σ**
  (×1.485 low; the face/side/back σ still uses it); **bin-midpoint face/side weights** (side
  over-weighted 6.8% in the outermost μ bins; replaced by exact bin averages).
- **GEANT4 option4 references** for this tier (GEANT4 relative to CeeLo, 15 m world): −4.87%
  (z −4.96) at 122 keV, of which 2.19 points are the Compton-model difference between the
  lists. **A 15 m GEANT4 world**: −3.02% (z −2.39) at 122 keV and −2.07% at 662 keV with
  Livermore. **Un-walked 511 photons**: +8.31% (z +6.31) at 2614 keV. **No bremsstrahlung
  channel**: a further +4.76% (z +2.20) there. Excluded by measurement while finding these:
  soil cross sections (all four components ≤0.16% vs a `G4EmCalculator` export, 20–300 keV
  — check μ_cs and μ_rs separately, a compensating pair is invisible in the total), the
  Compton angular distribution (CeeLo's tabulated S(x,Z) vs `G4LivermoreComptonModel`'s
  fitted form, 0.2% on soil; GEANT4 does not read the tabulated scattering functions), and
  the estimator's multiple-scatter accumulation (+0.15%, z +0.65, bin by bin against
  `run_brute`).
- **A room as a superposition of independently baked surfaces**: plus a ≥2-bounce diffuse
  bath +31%; plain superposition ±3% at 662 keV but +19.3% at 122 keV (an energy-dependent
  cancellation, not a model); one-bounce terms alone +22% to +62%. Infinite-plane surface
  kernels are the error. A scalar build-up multiplier has the wrong shape; 4V/S alone is not
  a room descriptor (two rooms at equal 4V/S differ by 10.9%). The finite-face slab-albedo
  operator (TODO.md) is the replacement design.
- **Tabulated (distance, height, E0) kernel grids with interpolation** (0.3%-precision
  nodes): in distance, log-linear interpolation misses by up to 12.6% and PCHIP on ln ŝ vs
  log d by 1.38%; in height the TOTAL interpolates to 0.55% while the 250–500 keV band, 39%
  of the flux, misses by up to 5.0% (z 8.2) — the total is a blind proxy for grid adequacy;
  across E0 band-by-band interpolation fails in absolute energy (−25.5%, z −38.7) and in
  E′/E0 (the backscatter edge slides across 1.6 bands). Replaced by the fit-first campaign
  on the registered energy axis.
- **Antipodal/hemisphere primary pairing and near-field cone antithetics**
  (`ExperimentalPrimarySampling`, `ExperimentalConeSampling`): do not meet the paired-bias
  and variance × wall-time gates; benchmark only.

## GEANT4 Validation

### Harness

GEANT4 11.4.0 is needed only for validation (source `bin/geant4.sh` before any run; CMake
config in `lib/cmake/Geant4`). The harness, `tools/geant4_validation/`, builds either
through `-DCEELO_BUILD_GEANT4_VALIDATION=ON` or standalone:

```bash
source <geant4-install>/bin/geant4.sh
cmake -B <g4val-build> -S tools/geant4_validation -DGeant4_DIR=<geant4-install>/lib/cmake/Geant4
cmake --build <g4val-build>
```

(On macOS, exit 134 / a dyld error means the GEANT4 libraries' Qt dependency is missing:
`brew install qtbase`.)

**Running a comparison:**
```bash
# 1. Export the GDML geometry + GPS macro from the engine, so GEANT4 runs exactly the MC's
#    geometry and source. benchmark_mc_configs exports configs 1-3, 5-7, 12 and 20-28
#    (`can_export`, a vetted opt-in list: an extended source's GPS volume must EXACTLY equal
#    the engine's source region); configs 8 and 11 are exported by generate_all_spectra.
cd build/examples
./benchmark_mc_configs --config 20 --energies 662,2614 --events 20000 --macro-events 12000000
#    -> detector_20.gdml + run_20_662keV.mac + run_20_2614keV.mac

# 2. Run GEANT4. Point sources: --cone-bias. Extended/shielded sources: isotropic, 16M+ events.
source <geant4-install>/bin/geant4.sh
<g4val-build>/ceelo_g4val detector_20.gdml run_20_662keV.mac g4val_cfg20_662keV.csv

# 3. Assemble tests/data/geant4_reference/<name>_multi.csv from the per-energy CSVs, add the
#    config to profiling/compare_validation.py CONFIGS, write our_N_multi.csv with
#    `benchmark_mc_configs --config N --precision p` (add --fep-window 0.75 to match a new
#    harness run), and compare:
python3 ../../profiling/compare_validation.py
```

**The exported source region is NESTED**, GEANT4's native idiom: each source volume is a
*full* solid at its own outer boundary carrying the next one in as a daughter, with only
the outermost placed in the world — no hollow solids, no boolean subtractions, no epsilon
between layers (two deliberate exceptions: the innermost sphere keeps `rmin = 1e-4` as the
GPS-at-the-centre workaround, and the Marinelli well-bottom annulus its 0.001 cm gap). A
point on a mother/daughter face is unambiguous to the navigator; a point on a shared face of
two siblings is not. It also keeps `/gps/pos/confine SrcMaterialPV` sampling the emitting
shell alone when the shell is cored (`IsSourceConfined()` locates the **deepest** volume).
Marinelli keeps its four-tubes-per-layer decomposition (an L-shaped fill is not concentric).

**Sphere angles are written at full precision**: `deltatheta` printed as `3.1415927`
exceeds π and makes a `G4Sphere` self-overlap (GEANT4 "pushes" stuck tracks). At 6M events
the pushes showed no measurable bias, so older references stand, but new ones must not be
generated on geometry that needs rescuing.

Output CSV columns: `energy_keV,fep_efficiency,fep_uncertainty,total_efficiency,total_uncertainty,num_events`.
Macros carry `/run/beamOn <N>` from `--macro-events`; `sed` it to change the count.

**Radioactive-ion (cascade) macros** (`/gps/particle ion` + `/gps/ion <Z> <A>`) — two
gotchas that silently give empty or wrong output:
- Add **`/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year`** after
  `/run/initialize`, or a long-lived parent (Ba-133, Cs-137, most calibration nuclides)
  never decays — the run "succeeds" with no counts.
- A **minutes-scale isomer is not fully decayed** under these settings (Cs-137 → Ba-137m:
  the 662 peak is ~2× low). Use a parent that decays directly (Ba-133).
- For electron-bearing channels (IC/β electrons crossing the gap) use an **air world**
  (`write_gdml(..., vacuum_world=false)`) so GEANT4's electrons air-degrade like the engine's.

**Harness flags** (`ceelo_g4val <gdml> <mac> <out.csv> [flags]`):
- `--cone-bias` — cone importance sampling (point sources only; matches the MC estimator).
  The cone is sized from the GDML's `ceelo_cone_bounds_cm` comment (older exports fall back
  to scanning the solids, which took z_min from the crystal face and so left an endcap or
  anything else in front of it partly outside the cone). Directions come from a generator
  seeded by `--seed`, so seeded processes are independent and their histograms can be
  summed.
- `--histogram` — write `<out>_histogram.csv` (`# Events:` header; lower-channel-energy
  bins; rows `energy_lower_keV,component,sum_w,sum_w2` at 12 significant digits — per
  emitted photon is `sum_w / Events`; bin 0 of `total` also holds the zero-deposit events);
  `--bin-width <keV>` sets the bin width (default 1 keV).
- `--correlated-gamma` — `G4PolarizationTransition` γ-γ angular correlation (cascade summing).
- `--diagnostics` — 11 per-interaction crystal-side histograms.
- `--entry-diag` — external-electron-entry counters at the crystal (u/s decomposition).
- `--lowcut` — production cut 0.001 mm (default 0.7 mm), for thin-detector secondaries.
  Fluorescence does not need it (deexcitation ignores the cut): on full spectra (option4,
  0.25 keV bins) bare NaI at 662 keV (4M) and NaI in a 1 cm Pb shell (12M) show no resolved
  difference with and without it (χ²/bin 1.03 and 1.13 on 4 keV bins, every region
  |z| ≤ 1.9, Pb K x-rays within 2%), at 6.4× and 2× the cost.
- `--em-standard` — `G4EmStandardPhysics`; `--em-option4` — `G4EmStandardPhysics_option4`.
- `--fep-precision <p>` — run until FEP relative precision ≤ p (cap 100M events).
- `--fep-window <keV>` — full-energy half-window (default `kDefaultFepWindowKeV`; 1.5 to match
  the older references).
- `--component-labels` — weighted histograms per interaction-history component (zero / FEP /
  escape / Compton / pair), for full-spectrum comparisons by family.
- `--sphere-flux R` — photon track-length flux φ = L/V in the scoring sphere, binned in
  (E, μ), all generations (the environment-tier observable).
- `--fluence` — fluence incident on the detector (primary only; air-scatter study).
- `--cone-scale <f>` — widen the auto cone half-angle by f (with `--cone-bias`); `--seed N` —
  the GEANT4 and cone-direction seed.
- `--vis`, `--vis-batch` — visualization.

Physics: the harness **default is `G4EmLivermorePhysics`** (`main.cc`); `--em-option4` selects
`G4EmStandardPhysics_option4`. **Every committed reference in `tests/data/geant4_reference/` was made
with option4, so pass `--em-option4` explicitly** when generating a new detector reference or
comparing against an existing one; otherwise the two runs use different EM models and the
discrepancy is an artifact. The environment-scatter tier is the exception: it is validated against
the Livermore default (see its section).

> **Reference physics vs CeeLo's data (open — TODO.md, "Generation/reference GEANT4 physics-list
> mismatch").** Two axes separate CeeLo from its detector references.
> - **Data evaluation:** CeeLo's photon magnitudes are EPICS2023; GEANT4 11.4.0 ships EPICS2017 (see
>   "EPICS2023 vs GEANT4's EPICS2017" below).
> - **Model set:** the references are option4. Option4 and Livermore share the photoelectric model
>   (`G4LivermorePhotoElectricModel`), the Rayleigh model, and the Klein-Nishina base for the Compton
>   total; they differ ONLY in the low-energy Compton model (`G4LowEPComptonModel` vs
>   `G4LivermoreComptonModel`). That difference is **angular, not attenuating**: 2.19% at 122 keV
>   and 0.39% at 2614 keV on a pure-scatter observable, while the *uncollided* term agrees to 0.44%
>   either way. Attenuation- and FEP-dominated results are largely insulated; **scatter-dominated ones
>   (total efficiency's scatter-in component, continuum shape, source-geometry leakage channels) are
>   the exposed set.** The 0.39% at 2614 keV is larger than the pooled ≥2 MeV FEP residual (+0.24%),
>   so only the experiment in TODO.md (one point-source config under both lists) bounds it.
> - Not every file in `tests/data/geant4_reference/` records its physics list in its header; new
>   references should.

Deexcitation: `SetFluo(true)`, `SetAuger(true)`, `SetDeexcitationIgnoreCut(true)` (the last is what keeps fluorescence
from heavy-element attenuators). To compare model choices read the GEANT4 source under
`source/processes/electromagnetic/` (`G4LowEPComptonModel`, `G4UrbanMscModel` /
`G4WentzelVIModel`, `G4SeltzerBergerModel`, `G4PEEffectFluoModel`, `G4Livermore*`), or add
`/run/verbose 1`, `/process/em/printParameters`, `/process/list` to a macro. The harness
records aggregate observables only; raw GEANT4 cross-section columns are not production
inputs and must not be committed.

### Validated Configurations

| Config | Detector | Attenuator | Source | FEP agreement | Total agreement |
|--------|----------|------------|--------|---------------|-----------------|
| 1 | 3"x3" NaI | bare | point, 10cm on-axis | ≤ 0.7% | ≤ 0.3% |
| 2 | 3"x3" NaI | 1mm Al | point, 10cm on-axis | ≤ 0.9%¶ | ≤ 0.4% |
| 3 | 2"x2" LaBr3 | 0.5mm Al | point, 5cm on-axis | ≤ 0.6% (≤ 1332 keV)¶; +0.8% @2000, +0.8% @3000 (LaBr3 residual, TODO.md)§ | ≤ 0.7%¶ |
| 5 | 1x1x0.5cm CZT | bare | point, 5cm on-axis | ≤ 0.6% (≤ 800 keV); −1.4% @1000, −2.2% @1500 (G4-reference-limited, \|z\| ≤ 1.6) | ≤ 0.6% |
| 6 | 3"x3" NaI | bare | point, 15cm 45deg off-axis | ≤ 0.5% | ≤ 0.4% |
| 7 | 3"x3" NaI | 1mm Al + 2mm Pb | point, 15cm on-axis | ≤ 0.7% (≥ 200 keV) | ≤ 0.4% (≥ 200 keV)*** |
| 8 | 3"x3" NaI | 0.5mm Al | Marinelli beaker, water | ≤ 0.8% (≥ 100 keV); −1.9% @ 59* | ≤ 0.4%* |
| 11 | 3"x3" NaI | bare | point, 10cm + 0.5cm Fe shield | ≤ 0.8% (200-3000)** | ≤ 0.7%** |
| 12 | 3"x3" NaI | bare | 10x15x20cm SS304 box, cellulose | ≤ 0.2% (≥ 200 keV); −3.9% @ 59 | ≤ 1.2% (≥ 200 keV); −2.9% @ 59 |
| 20 | 3"x3" NaI | bare | 10x15x20cm Pb-wall box, cellulose | ≤ 1.1% | ≤ 1.1%**** |
| 21 | 3"x3" NaI | bare | 10x15x20cm W-wall box, cellulose | ≤ 1.7% | +2.0% @662, +2.9% @2614 (TODO.md)**** |
| 25 | GEM35-70 HPGe coax, sharp edge | bare | point, 5cm on-axis | ≤ 0.9% | ≤ 0.2% |
| 26 | GEM35-70 HPGe coax, **bulletized** + round-tipped bore | bare | point, 5cm on-axis | ≤ 0.7% | ≤ 0.3%† |
| 27 | GEM35-70 HPGe coax, bulletized | bare | point, 2cm + 0.5cm Fe shell | +0.1% @122; −8.2% @60‡ | +0.2% @122‡ |
| 28 | GEM35-70 HPGe coax, bulletized | bare | point, 10cm + 0.5cm Fe shell | +0.4% @122; −8.7% @60‡ | −0.04% @122‡ |

CeeLo at `--precision 0.003` (`benchmark_mc_configs --config N`, all 126 config/energy
points; this branch's engine, regenerated Sep 2026; cfg 27/28 with `--fep-window 0.75`,
cfg 20/21 with `--energies 662,2614`) against the committed references, on the EPICS2023
photon data and after the crystal-walk fixes (Physics above).
A precision target is met only when *every* requested target is (FEP as well as total);
a few tenths of a percent of every entry is still statistics. `compare_validation.py`
prints every per-energy deviation with both 1σ and z. GEANT4 references: 64M events per
energy for cone-biased point sources, 16M isotropic for extended/shielded sources (64M for
cfg 11 at 2000/3000 and cfg 12 at 662/1173/2614 keV).

¶ **Two entries are set by single statistical draws**, confirmed by re-running just those rows at
0.1% precision: cfg 2 FEP @800 keV reads −0.85% (z −3.2) in the 0.3% sweep and +0.06% (z +0.5)
re-run; cfg 3 @30 keV reads FEP +0.50% / total +0.69% (z 3.1 / 4.6) in the sweep and +0.11%
(z +1.2) / +0.23% (z +2.6) re-run. Without those two rows cfg 2 FEP is ≤ 0.6% and cfg 3 total
≤ 0.3%. The table keeps the sweep's values; the harness is not bit-reproducible run to run.

§ **The LaBr3 residual is confined to ≥ 2 MeV.** At 1000 / 1173 / 1332 keV config 3 runs
−0.23 / +0.17 / +0.60% (|z| ≤ 2.0). What survives is the ≥2 MeV structure the TODO item describes:
**+0.81% (z 2.6) @2000 and +0.77% (z 2.4) @3000**. Those z-scores are *lower* than the +4.0 / +2.6
quoted in TODO.md only because that item was measured in a dedicated 0.16%-precision study, not at
the 0.3% these rows carry — the deltas agree, so the item stands.

**Configs 25/26 are a matched pair**: the same GEM35-70 coax with a sharp and a
bulletized front edge. The quantity the pair measures is the **bulletization effect**,
Δ = 1 − ε(26)/ε(25), computed independently in each code (values as
`compare_validation.py` prints them):

| E (keV) | Δ FEP CeeLo | Δ FEP G4 | difference | z |
|---------|-------------|----------|------------|---|
| 45 | 16.37% | 16.39% | −0.02 pp | −0.10 |
| 59.5 | 15.20% | 15.19% | +0.01 pp | 0.07 |
| 88 | 12.21% | 12.25% | −0.04 pp | −0.20 |
| 122 | 8.72% | 9.09% | −0.36 pp | −1.68 |
| 344 | 3.84% | 4.09% | −0.25 pp | −0.70 |
| 662 | 2.97% | 3.15% | −0.18 pp | −0.45 |
| 1332 | 2.70% | 2.47% | +0.23 pp | 0.53 |

Every point agrees within |z| ≤ 1.7 on an effect running from 16% down to 2%. The effect is
far larger than the 1.28% of crystal volume the fillet removes: from this source position
bulletization removes 17% of the rays that hit the crystal but only 2.4% of the path
length (it eliminates short grazing corner clips), which costs most at low energy where
photons stop in the first millimetres. A dedicated ~0.03%-precision study found a small
shared germanium offset (+0.135% at 88 keV in *both* configs, ~−0.1% at 122 keV) that
cancels in Δ (TODO.md). The pair validates "CeeLo's tracer agrees with CeeLo's exported
solid under GEANT4 physics"; the link to ANGLE's `bulletizingRadius` rests on the Pappus
closed form and the unit tests in `tests/test_geometry.cpp` and
`tests/test_geant4_export.cpp`.

‡ **Configs 27/28 are the deep-iron low-energy probe** (GEM35-70 + point source + 0.5 cm
Fe shell at 2 / 10 cm). The absolute 60 keV FEP gap is the EPICS2023-vs-EPICS2017 iron
photoelectric difference (~2% in μ, CeeLo = NIST) through the shell's 4.6 mfp, the same
at both distances as a cross-section effect must be; 60/88 keV are `SKIP`ped and 122 keV
gates at 1%. The u/s decomposition on both sides agrees on the *scattered share* of the
peak (60 keV: 30.0% vs 30.6% at 2 cm), so the forward-Rayleigh / small-angle-Compton
placement behind thick iron matches GEANT4 to ~2% of that stream.

\* **Config 8 (water Marinelli)**: the reference is a 32M-event run on the current
geometry, gating all six energies (FEP tolerance 2.5%, total 1.0%). The residual at
59 keV is the scattered FEP stream: the unscattered FEP `eps_u` matches GEANT4 to
0.01–0.3% at every energy, but the small scattered stream `eps_s` (forward Rayleigh /
small-angle Compton landing in the window) is ~15% low. Its electron escape uses the
regime-aware `walk_in_source_geometry` (water → fully analog).

\*\* **Configs 11/12**: the cfg 11 GDML uses `rmin = 1e-4 cm` for the Fe shell (GPS at the
exact centre of an `rmin = 0` sphere misnavigates). Total agreement above 1022 keV requires
the source-shield secondary channels (Source shielding above), up to ~6% of total at
3 MeV; the electron-escape channel must stay on (without it cfg 11 @3000 total is −3.05%,
z 17). The escape is calibrated by an analog surface-emergence (1 − albedo) gate on the
escape outcome (brems kept): it suppresses the sub-MeV skin escapers GEANT4 backscatters
and saturates to ≈1 above ~1.5 MeV (Known differences → source-electron skin escape).

\*\*\* **Config 7 at 100 keV** (~−25 to −31% total) is excluded from the gate (`SKIP`).
GEANT4's *in-tracking* μ_total(Pb) is smeared across the 88 keV K edge relative to its own
`G4EmCalculator` export (which matches XCOM to 0.23%): measured through the 2 mm Pb stack,
Δμ/μ = +0.35% at 80 keV, −18% at 90 keV, −4.4% at 100 keV — the classic λ-table edge
interpolation signature. CeeLo's tables equal XCOM and below-edge transmission matches
GEANT4. (A follow-up — CeeLo's escaped Pb Kα/Kβ ratio 0.67 vs GEANT4's 0.82 — is in TODO.md.)

\*\*\*\* Configs 20/21 (Pb/W-wall boxes) have 16M isotropic option4 references at 662 and
2614 keV only, ~1% relative at 2614 keV, so a ~1% CeeLo shift there sits inside the gate's 3σ
allowance and the gate alone cannot see it. They do not show a common high-Z total offset at
2614 keV: a hypothesised +2–3.7% excess in both escape modes does not reproduce (Pb −0.10% ± 0.53%,
z −0.19; W +1.45% ± 0.56%, z +2.58 before the crystal-walk changes, ≤ 0.9% in the table above).
Pb's residual is at 662 keV. Emitting the walls' K x-rays and transporting source
positrons adds ~+1% of total at 662 keV and +0.5–0.7% at 2614 keV to both boxes
(A/B with `--no-source-fluorescence` / `--no-source-positrons`): Pb's 662 keV total moves
from −0.96% to +0.26%, and W's existing excess grows to +2.0% / +2.9% (still inside the
gate's tolerance + 3σ).

**The CeeLo side of the gate is whatever `our_*_multi.csv` sits in `build/examples/`.**
`benchmark_mc_configs --config N --precision p` writes those files locally (they are not committed;
`tests/data/ceelo_reference/README.md`) and `profiling/compare_validation.py` reads them from there
(`--mc-dir` overrides). They must be regenerated whenever the engine changes, or the gate reports
staleness as a physics discrepancy — the cfg-8 low-energy deficit and a "≥2 MeV FEP family"
(mean −0.78%, −0.06% once regenerated) were stale-reference artifacts, while the Na-22 P511
drift was real and went unnoticed (TODO.md). Nothing
checks their age and the gate is **not** wired into ctest — see `TODO.md`.

**EPICS2023 vs GEANT4's EPICS2017.** CeeLo's photon magnitudes come from a *newer*
evaluation than GEANT4 11.4.0 ships, so exact agreement is not the expected outcome.
Migrating from the previous GEANT4-fitted representation moved the mean |deviation| over
174 gated points from 0.366% to 0.432% (70 closer, the rest further); the largest moves were
config 12 @59 keV (see Known differences) and config 11 @3000 keV. Acceptance: 0.72% p95 /
1.64% worst contribution-weighted attenuation change, with independent NIST XCOM checks
favouring the direct EPICS values
(`tools/prepare_cross_sections/reports/migration_acceptance.json`).

### Spherical / annular source geometries

`examples/sphere_g4_validation.cpp`: three 3"×3" NaI geometries — **G-A** bare
self-attenuating thorium sphere (R = 2 cm), **G-B** trace source in a soil sphere
(R = 3 cm) + Al(0.2)/Fe(0.5) shells, **G-C** void-centre soil shell [2,3] cm + Fe(0.5) — at
~2 cm and ~50 cm, Th-232 chain energies, 4M-event isotropic GEANT4 volume sources.

- **FEP agrees everywhere**: |z| ≤ 1.4 at 238.6/583/911 keV, validating geometry,
  sampling and self-attenuation.
- **Low-Z totals agree**: G-C, re-measured with source-electron transport on and a
  0.75 keV window, agrees to ≤0.74% (source-core table below). (G-B was measured with the
  channel off, which costs several percent of total above 1022 keV.)
- **High-Z self-attenuating total**: G-A (near, against 16M-event option4 runs) reads total
  −0.12 ± 0.23% / +0.03 ± 0.19% / +0.10 ± 0.16% / +0.01 ± 0.10% (z −0.5 / +0.2 / +0.6 /
  +0.1) at 238.6/583/911/2614.5 keV and FEP −0.22 / +0.05 / +0.63 / +0.14 ± 0.28% (z −0.8 /
  +0.2 / +2.2 / +0.5), with source K and L x-rays. Before source K fluorescence was emitted
  the totals read −22% / −4.1% / −2.6%. Th L x-rays (13–19 keV) leave only from the
  sphere's skin: emitting them moves the totals by under 0.06%.

At ~50 cm the MC converges well via the two-stream direct stream; a matching GEANT4 run
needs much higher isotropic statistics and has not been done. The cascade path on a
sphere is exercised by `tests/test_spherical_source.cpp`.

### Source cores — attenuating layers inside a hollow source

Three geometries at the NEAR distance, all with an **iron** core (its K x-rays reabsorb
within microns, keeping TOTAL judgeable), plus G-C as the uncored control: **G-D** (G-C's
shell with an Fe-filled centre), **G-F** (closed-cavity cylinder, R=3/hz=3, cavity r=2/hz=2,
soil + Fe core), **G-G** (box shell, half-dims 3 cm, cavity 2 cm, soil + Fe core). 16M
isotropic GEANT4 events per energy, CeeLo at 0.2% on FEP, both at `kDefaultFepWindowKeV`,
vacuum world, source-electron transport on.

| Config | FEP agreement | Total agreement |
|--------|---------------|-----------------|
| G-C (uncored control) | −0.55 to +0.08% (\|z\| ≤ 2.3) | −0.63 to +0.74% (\|z\| ≤ 3.6) |
| G-D cored sphere | −0.30 to +0.19% (\|z\| ≤ 1.0) | −0.00 to +0.39% (\|z\| ≤ 1.6) |
| G-F cored cylinder | −0.12 to +0.09% (\|z\| ≤ 0.5) | −0.14 to +0.47% (\|z\| ≤ 2.4) |
| G-G cored box | −0.74 to +0.03% (\|z\| ≤ 2.7) | −0.36 to +0.26% (\|z\| ≤ 1.8) |

Largest deviation over the 16 (case, energy) rows: 0.74%, cored and uncored alike. FEP is
statistically clean (χ²/n = 1.1); TOTAL is not (χ²/n ≈ 3.2, positive mean, same in cored
and uncored cases — a property of the extended-source path, tracked in TODO.md). The iron
core removes 10–18% of the peak (238.6 → 2614.5 keV), so a missing or mis-sized core would
be obvious.

Two export controls exist because the GEANT4 geometry is generated by CeeLo:
- **G-E, self-core vs solid** (a shell cored with its own material is optically a solid
  ball but emits only from the shell): FEP ratio CeeLo 1.0330 ± 0.0028 vs GEANT4
  1.0354 ± 0.0022 (z −0.66), the GEANT4 ratio 16σ from the 1.0000 that GPS leaking into
  the core would give. It cannot detect a mis-sized core (the core is optically
  invisible); G-D/F/G and the GDML-vs-tracer tests cover that.
- **G-Eadd, additivity** (one 2 cm core vs four 0.5 cm cores): CeeLo 1.0045 ± 0.0028
  (+1.6σ). Only CeeLo's arm is informative (GEANT4 sees the same soil either way).

The nested export itself showed no measurable physics change against the pre-nesting GDML
(G-C at 2614.5 keV, 16M each: FEP z −1.21, total z −0.30 — a bound of ±0.9% / ±0.3%);
sub-percent assurance comes from the GDML-vs-tracer tests, which compare the exported
solids against `trace_source_segments` exactly.

### Generating Reference Data

```bash
cd build/examples
./generate_all_spectra [n_mc_histo] [n_g4_macro] [n_mc_multi]
# Defaults: 4M histogram, 100M G4 macro, 100k efficiency; scores at kDefaultFepWindowKeV
bash run_g4_all.sh    # generated script handles cone bias per config
```

`ceelo_g4val ... --diagnostics` gives per-interaction GEANT4 histograms for debugging
scatter angles, interaction depths and escape energies.

## Histogram CSV Format

```
# MC energy deposit histogram
# Config description
# Primary energy: 662 keV
# Events: 4000000
# Bin width: 1 keV
RealTime: 123.456 s
Energy (keV),Counts
0.0,1234
1.0,5678
...
```

- Lower channel energy (not bin centre).
- The `RealTime` line has no `#` prefix (spectral analysis software reads it):
  `N_events / (1E6 * omega_frac)` cone-biased, `N_events / 1E6` isotropic — the equivalent
  measurement time for a 1 MBq source with 100% branching ratio.

## Known differences and approximations

Standing explanations that are easy to trip over. Open work on any of them is in TODO.md.

- **Config 12 at 59 keV is a data difference, not a transport defect.** Its FEP moved
  from −0.05% to −3.76% (z 8.9) across the EPICS2023 migration. μ at 59 keV changed by
  +1.73% for SS304, +0.25% cellulose, −0.01% NaI; the 0.2 cm steel wall has τ = 1.99, so
  first-order attenuation predicts −3.79% against −3.72% ± 0.42% observed, and the same
  model reproduces the energy dependence (≈0 above 200 keV). The change is a smooth
  ~30 keV-wide bump in iron's photoelectric cross section peaking at +2.23% near 58 keV.
  Arbitrated against NIST (Hubbell & Seltzer SRD 126, derived from neither evaluation):

  | element | 50 keV | 60 keV | 80 keV | 100 keV |
  |---|---|---|---|---|
  | Fe pre-migration | −1.30% | −1.26% | −0.36% | +0.08% |
  | Fe EPICS2023 | **+0.31%** | **+0.36%** | **+0.39%** | **+0.32%** |
  | Cr pre-migration | −1.67% | −1.32% | −0.34% | +0.11% |
  | Cr EPICS2023 | **+0.31%** | **+0.35%** | **+0.33%** | **+0.27%** |
  | Ni pre-migration | −0.97% | −1.07% | −0.42% | +0.06% |
  | Ni EPICS2023 | **+0.31%** | **+0.38%** | **+0.35%** | **+0.29%** |

  EPICS2023 is ~4× closer to NIST at 50–60 keV and its offset is flat, the signature of a
  coherent evaluation. **So config 12 @59 keV is CeeLo being more accurate than GEANT4.**

- **Config 8 at 59 keV is separate**: water's μ moved only +0.26% and the beaker's τ ≲ 1,
  so attenuation explains ≲0.3% of it; it is the scattered-FEP stream (footnote *).

- **Doppler broadening** (impulse approximation, default ON) broadens the Compton *edge*
  to match GEANT4; the subshell binding energy is deposited locally (or, with spectrum
  fidelity, a K vacancy relaxes), so FEP and total are statistically unchanged. **Never cite Doppler broadening as a cause of FEP/total
  discrepancies.**

- **Source-electron skin escape** uses an empirical surface-emergence albedo gate
  (`ElectronCsda::source_escape_survival_exit`), evaluated at the EXIT state (exit KE +
  exit-material Z) with a light-Z floor (η₀ ≤ 0.16 → fully analog, so water and Al
  positrons match GEANT4 with no per-caller flag) and an exit-energy window. Validated
  across Z = 5–82 (PE/water/Al/Fe/SS304/Sn/W/Pb) to ~1–2%. The model is **compile-time
  selectable** (`-DCEELO_SOURCE_ESCAPE_MODEL=gate|tails|gs`, default `gate`; env A/B knobs
  `MCDET_NO_{STRAGGLE,B2,TAILS}`): `tails` / `gs` are first-principles (straggling +
  screened-Rutherford hard tail, no gate), differing only in the soft Gaussian-core width,
  and only inside `walk_in_source_geometry`.

  | Model | Gaussian core θ₀ | hard tail | straggling | escape decision |
  |-------|------------------|-----------|------------|-----------------|
  | **gate** | Highland (per-step) | — | — | empirical exit-state (1−albedo) gate |
  | **tails** | Highland − hard variance | screened-Rutherford μ<μ_c | Bohr/clamped | analog |
  | **gs** | √(G₁,soft·Δx) (screened-Rutherford soft moment) | screened-Rutherford μ<μ_c | Bohr/clamped | analog |

  Total-efficiency error vs GEANT4 at 662 keV (%, MC ~0.2%):

  | Config | Z / geom | gate | tails | gs |
  |--------|----------|------|-------|-----|
  | cfg 24 PE    | ~5  box   | −0.80 | −0.57 | (≈−0.6) |
  | cfg 8 water  | ~7.5 L-shape | −0.95 | (≈−1.0) | −1.06 |
  | cfg 23 Al    | 13  box   | −0.39 | −0.36 | −0.57 |
  | cfg 11 Fe    | 26  shell | −0.17 | +0.13 | +0.13 |
  | cfg 12 SS304 | 26  box   | **+0.46** | +1.07 | +1.65 |
  | cfg 22 Sn    | 50  box   | **+0.32** | +1.17 | +1.69 |
  | cfg 21 W     | 74  box   | +0.95 | +1.94 | (high) |
  | cfg 20 Pb    | 82  box   | −1.46 | **−0.36** | −0.58 |

  Mean |err|: gate 0.56% < tails 0.62% < gs 0.92%. `tails` is better at extreme Z but
  over-diffuses sub-MeV mid-Z **boxes** (a large-angle scatter is as likely to redirect
  toward a wall as away); on spheres/shells all three agree to ~±1%. `gs` under-scatters
  ~30% (its screened-Rutherford G₁ lacks the Mott and nuclear-form-factor corrections
  GEANT4 uses). High Z is not one knob: Sn/W over-predict while Pb under-predicts at the
  same capped albedo, and Pb's −1.46% is mostly photon transport (Pb stops electrons). Run
  time is immaterial (`tails` 1.03–1.10×, `gs` 1.1–1.2× `gate`).

- **Electron tables above uranium are uranium's.** ESTAR stopping powers, ICRU-49 mean
  excitation energies and the EPQ Seltzer–Berger tables stop at Z = 92, so Np–Cf are
  treated as uranium throughout the electron physics (`electron_table_z()`): the tables
  and the atomic weight, Z/A, Z(Z+1)/A, radiation length and effective Z used with them,
  keeping each material's electron view self-consistent. The photon side uses each
  element's own data and mass. Against ESTAR (per unit mass; negative = reusing U
  underestimates):

  | Z | I (eV; U 890) | collision | radiative | CSDA range |
  |---|---:|---:|---:|---:|
  | 93 Np | 902 | −1.2 .. −1.0% | −2.9 .. −1.9% | +0.3 .. +1.6% |
  | 94 Pu | 921 | −1.3 .. −0.4% | −4.5 .. −2.6% | −1.3 .. +2.1% |
  | 95 Am | 934 | −1.8 .. +0.8% | −5.2 .. −2.4% | −3.3 .. +2.1% |
  | 96 Cm | 939 | −1.1 .. +1.6% | −6.0 .. −2.2% | −4.3 .. +1.7% |
  | 97 Bk | 952 | −1.9 .. +1.1% | −8.2 .. −3.6% | −4.7 .. +2.9% |
  | 98 Cf | 966 | −2.1 .. +2.3% | −8.9 .. −3.4% | −6.8 .. +2.8% |

  The largest range errors are at the 10–20 keV nodes (micron ranges). The efficiency
  effect has not been measured; the reuse is pinned in `test_csda`
  (`ActinideElectronReuse`).

- **Above 3 MeV transport is not GEANT4-validated.** The photon tables reach 20 MeV so
  attenuation is available for high-energy reaction lines; a narrow-beam PuO₂/Pb test at
  15 MeV agrees with exp(−μL) to 0.09%. Pair production splits the kinetic energy evenly
  between the leptons, and photonuclear absorption is not modelled (as in EPDL and XCOM;
  a few percent of Pb/U attenuation at the 12–15 MeV giant resonance).

- **HPGe at 59.5 keV behind a thick dead layer is a data difference too.** On
  `detective_x_fit` (0.0739 cm front dead layer, τ ≈ 0.82) CeeLo reads ~1.0–1.1% below
  GEANT4 in FEP and total at 59.5 keV with either GEANT4 list. GEANT4's photoelectric
  cross section for light and mid-Z elements carries the same smooth deficit near 60 keV
  found for iron above: CeeLo / GEANT4 at 59.54 keV is Ge +1.25%, Al +2.1%, C/O +2.2%
  (I, Cd, Te, Pb within 0.15%; `studies/xs_probe/ge_59keV.md`). Against NIST (SRD 126),
  μ/ρ(Ge) is CeeLo +0.22…+0.43% at 30–100 keV (flat) and GEANT4 +0.9% at 30–40 keV,
  −0.72% at 60 keV: CeeLo is the closer. Through the dead layer and endcap that predicts
  −0.97% in transmission.

- **Pair production near threshold is a GEANT4 option4 difference.** option4's
  `G4BetheHeitler5DModel` takes its cross section from `G4PairProductionRelModel::
  ComputeParametrizedXSectionPerAtom`, a fit valid from 1.5 MeV that below it scales the
  1.5 MeV value by ((E − 2mc²)/(1.5 MeV − 2mc²))². Against NIST XCOM (nuclear + electron
  field, the local XCOM database; `studies/xs_extend/pair_threshold_xcom.py`), EPICS2023
  agrees to 4 digits at every XCOM node from 1.25 to 3 MeV for Na, Al, Ge, Cd, Te, I and
  Pb, and GEANT4 Livermore (`G4EMLOW8.8/epics2017/pair`) equals EPICS2023 to <0.1%. option4
  is the outlier:

  | option4 / XCOM | 1250 keV | 1500 keV | 2000 keV | 3000 keV |
  |---|---|---|---|---|
  | Na | 1.199 | 0.959 | 0.955 | 1.004 |
  | Ge | 1.212 | 1.075 | 1.020 | 1.006 |
  | I | 1.121 | 1.065 | 1.026 | 1.005 |
  | Pb | 1.045 | 0.962 | 0.963 | 0.986 |

  (At 1332 keV option4 is +3–10% above EPICS2023 for Z = 11–53, Ge +9.9%, and −4.9% for Pb.) The single
  and double escape peaks at 1332 keV therefore read 5–13% low against option4 but agree
  with Livermore: HPGe SE 1.020 ± 0.020 (z +1.0), DE 0.990 ± 0.025 (z −0.4) (full-spectrum
  study). Between XCOM's sparse nodes its own spline interpolation differs from EPICS2023's
  dense table by −7% to +3% (1.1–1.35 MeV); quote XCOM only at its nodes.

## Common Pitfalls

- **Two Z domains**: `kMaxZ` (98) is the photon and material domain; electron tables stop
  at `kMaxElectronTableZ` (92). Index an electron table only through `electron_table_z()`.
  The hand-maintained `photon_epics_data.h` and `element_data.h` are sized by these
  constants and the generated sources `static_assert` against them, so a changed domain
  fails to compile until the tables are regenerated.
- **Physics tables are generated**: do not edit `element_data.cpp`, `photon_epics_data.cpp`,
  `relaxation_epics_data.cpp` or `estar_stopping_data.cpp` by hand. Use the locked
  generators and update the generated manifest last.
- **C++ `const` linkage**: namespace-scope `const` arrays have internal linkage; if a
  header declares `extern const T arr[]`, the definition must also say `extern`.
- **Eigen `.cross()`** needs `#include <Eigen/Geometry>`, not just `<Eigen/Core>`.
- **Batch size**: the convenience `compute(double, …)` uses `max(10000, num_events/10)` so
  threads check termination often; `compute(SimulationConfig)` uses `config.batch_size`
  (default 10000, floor 100) for the first wave and doubles later waves.
- **Importance-sampling denominator**: `(1/N) * sum(w_i * I_i)`. Using `sum(w_i)` gives
  intrinsic efficiency (it cancels the cone weight), not absolute efficiency.
- **Code placement moves timings by ~10%**: the same library code linked at different
  addresses ran `benchmark_mc_configs` 8–11% slower in CPU time on Apple silicon whenever
  the Fast electron walk landed differently (a driver that only parsed three new flags did
  it). Before attributing a cost change to a code change, rebuild the comparison from the
  same driver and check the hot functions' addresses (`nm`); measure CPU time, interleaved,
  at fixed event counts.
- **Coincident GEANT4 boolean surfaces** stall navigation ("track stuck") and — worse —
  inflate the `--entry-diag` electron-entry rate ~1.75× while leaving totals unchanged, so
  they silently corrupt channel-level comparisons. The nested source export avoids them;
  keep it that way. (~1.5k benign warnings per run at `active_crystal_PV` remain in every
  config.)
