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

The goal is agreement with the large, fully detailed Monte Carlo engines (in practice
GEANT4 11.4.0) to ~1% from 30 keV to 3 MeV, with anything beyond that understood and
documented (Validated Configurations and Known differences below). Full-energy peak and total
efficiencies are the priority; capturing the spectrum accurately is also wanted,
and is worth a few percent of slowdown but not much more.

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
  geometry/         Ray tracing (cylinder, box, sphere), attenuator shells, source geometry
  io/               Storable detector-response payload, MC-free response kernel,
                    EFFTRAN-style efficiency transfer (CEELO_BUILD_RESPONSE)
  materials/        Material definitions (composition, density), builtin materials
  physics/          Electron CSDA range tables, stopping powers, Moliere scattering, FEP window
  transport/        Photon transport (Compton, PE, Rayleigh, PP, fluorescence, bremsstrahlung)
examples/           Examples; benchmark_mc_configs + bench_configs.h (validation configs)
tests/              Boost.Test suites + reference data
tools/
  prepare_cross_sections/  Locked EPQ/EPICS2023/ESTAR generators and reports
  geant4_validation/       GEANT4 harness (ceelo_g4val); generate_all_spectra (built with examples)
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
  (radioactive-decay daughters one beyond the photon tables). L lines are deposited
  locally in the crystal (TODO.md). The same K relaxation (`sample_k_relaxation`)
  follows photoelectric absorption in **source material and shields**
  (`set_source_fluorescence`, default on; the x-ray joins the source secondaries) and,
  with spectrum fidelity, the K vacancy a Doppler-sampled **Compton** collision leaves
  (`TransportConfig::enable_compton_fluorescence`).

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
that shapes the pulse-height spectrum but moves FEP and total by well under 0.1%, so
FEP/total runs leave it off (the default). Here it sets
`TransportConfig::enable_compton_fluorescence` (`enable_compton_fluorescence`): the K
vacancy a Doppler-sampled Compton collision leaves (`ComptonScatterResult::shell_index`)
relaxes through `sample_k_relaxation` and its x-ray is transported as after photoelectric
absorption (`FluorescenceRecord::from_compton`); off, the whole binding energy is
deposited on the spot. It restores the K-escape peaks above ~300 keV (NaI, CZT: 0.90–0.95
of GEANT4 → 1.00–1.02).

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
  source material/shields; 511 keV pairs from pair production there; bremsstrahlung from electrons stopping there (Seltzer–Berger along a
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
  once cost −26% at 30 keV). A crystal with a **bore and/or bulletized edge** is one
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
0.75 keV (−3.09%, z −13.5); bare-NaI config 1 shows none (|z| ≤ 1.0). The harness
(`EventAction`) and `generate_all_spectra` already score at 0.75 keV, with no CLI flag to
change the harness, so **their new output is not comparable with the 1.5 keV references**.

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
  fixed-grid only; `closed_loop` refuses it. There is no near-field table behind the face plane:
  the kernel carries the distance dependence there.
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
  points. `back_structure_approx` (0.10) is the cost of a reasonable back description only:
  leaving 2 cm of steel out of an HPGe's back changed its efficiency directly behind 30-fold at
  122 keV.
- **Serialization:** a response with back nodes is stamped version 4, which older readers refuse
  (they would interpolate across the seam); everything else keeps its version and hash.
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
Rayleigh); it is not yet validated against GEANT4, and air scatter is not modelled.

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
- `--histogram` — write `<out>_histogram.csv` (`# Events:` header; lower-channel-energy
  bins); `--bin-width <keV>` sets the bin width (default 1 keV).
- `--correlated-gamma` — `G4PolarizationTransition` γ-γ angular correlation (cascade summing).
- `--diagnostics` — 11 per-interaction crystal-side histograms.
- `--entry-diag` — external-electron-entry counters at the crystal (u/s decomposition).
- `--lowcut` — production cut 0.001 mm (default 0.7 mm), for thin-detector secondaries;
  fluorescence is already kept by `SetDeexcitationIgnoreCut(true)`.
- `--em-standard` — `G4EmStandardPhysics` instead of option4.
- `--fep-precision <p>` — run until FEP relative precision ≤ p (cap 100M events).
- `--vis`, `--vis-batch` — visualization.

Physics: **FTFP_BERT + `G4EmStandardPhysics_option4`**; deexcitation `SetFluo(true)`,
`SetAuger(true)`, `SetDeexcitationIgnoreCut(true)` (the last is what keeps fluorescence
from heavy-element attenuators). To compare model choices read the GEANT4 source under
`source/processes/electromagnetic/` (`G4LowEPComptonModel`, `G4UrbanMscModel` /
`G4WentzelVIModel`, `G4SeltzerBergerModel`, `G4PEEffectFluoModel`, `G4Livermore*`), or add
`/run/verbose 1`, `/process/em/printParameters`, `/process/list` to a macro. The harness
records aggregate observables only; raw GEANT4 cross-section columns are not production
inputs and must not be committed.

### Validated Configurations

| Config | Detector | Attenuator | Source | FEP agreement | Total agreement |
|--------|----------|------------|--------|---------------|-----------------|
| 1 | 3"x3" NaI | bare | point, 10cm on-axis | ≤ 0.6% | ≤ 0.5% |
| 2 | 3"x3" NaI | 1mm Al | point, 10cm on-axis | ≤ 0.4% | ≤ 0.5% |
| 3 | 2"x2" LaBr3 | 0.5mm Al | point, 5cm on-axis | ≤ 0.8%; +0.7% @2000, +0.6% @3000 (TODO.md) | ≤ 0.4% |
| 5 | 1x1x0.5cm CZT | bare | point, 5cm on-axis | ≤ 0.8% (≤ 800 keV); −1.5% @1000, −2.5% @1500 (G4-reference-limited, \|z\| ≤ 1.7) | ≤ 0.5% |
| 6 | 3"x3" NaI | bare | point, 15cm 45deg off-axis | ≤ 0.5% | ≤ 0.8% |
| 7 | 3"x3" NaI | 1mm Al + 2mm Pb | point, 15cm on-axis | ≤ 0.8% (≥ 200 keV) | ≤ 0.5% (≥ 200 keV)*** |
| 8 | 3"x3" NaI | 0.5mm Al | Marinelli beaker, water | ≤ 0.9% (≥ 100 keV); −1.6% @ 59* | ≤ 0.6% |
| 11 | 3"x3" NaI | bare | point, 10cm + 0.5cm Fe shield | ≤ 0.5% (200-3000)** | ≤ 0.8%** |
| 12 | 3"x3" NaI | bare | 10x15x20cm SS304 box, cellulose | ≤ 0.2% (≥ 200 keV); −4.0% @ 59 | ≤ 1.4% (≥ 200 keV); −3.2% @ 59 |
| 25 | GEM35-70 HPGe coax, sharp edge | bare | point, 5cm on-axis | ≤ 0.5% | ≤ 0.4% |
| 26 | GEM35-70 HPGe coax, **bulletized** + round-tipped bore | bare | point, 5cm on-axis | ≤ 1.1% | ≤ 0.4% |
| 27 | GEM35-70 HPGe coax, bulletized | bare | point, 2cm + 0.5cm Fe shell | +0.1% @122; −8.2% @60‡ | +0.2% @122 |
| 28 | GEM35-70 HPGe coax, bulletized | bare | point, 10cm + 0.5cm Fe shell | +0.4% @122; −8.7% @60‡ | −0.0% @122 |

CeeLo at `--precision 0.003` (`benchmark_mc_configs --config N`; Sep 2026, main @21d883e)
against the committed references, on the EPICS2023 photon data and after the crystal-walk
fixes (Physics above). Later main changes (photon tables extended to Cf, the source-core
fixes) move μ by ≤0.1% and do not touch these geometries. Configs 27/28 are scored at
0.75 keV on both sides (`--fep-window 0.75`), from a Sep 22 regeneration.
A precision target is met only when *every* requested target is (FEP as well as total);
a few tenths of a percent of every entry is still statistics. `compare_validation.py`
prints every per-energy deviation with both 1σ and z. GEANT4 references: 64M events per
energy for cone-biased point sources, 16M isotropic for extended/shielded sources (64M for
cfg 11 at 2000/3000 and cfg 12 at 662/1173/2614 keV).

**Configs 25/26 are a matched pair**: the same GEM35-70 coax with a sharp and a
bulletized front edge. The quantity the pair measures is the **bulletization effect**,
Δ = 1 − ε(26)/ε(25), computed independently in each code (values as
`compare_validation.py` prints them):

| E (keV) | Δ FEP CeeLo | Δ FEP G4 | difference | z |
|---------|-------------|----------|------------|---|
| 45 | 16.34% | 16.39% | −0.05 pp | −0.32 |
| 59.5 | 15.50% | 15.19% | +0.32 pp | 2.06 |
| 88 | 12.22% | 12.25% | −0.03 pp | −0.16 |
| 122 | 9.27% | 9.09% | +0.19 pp | 0.88 |
| 344 | 4.51% | 4.09% | +0.42 pp | 1.20 |
| 662 | 3.01% | 3.15% | −0.14 pp | −0.36 |
| 1332 | 1.82% | 2.47% | −0.65 pp | −1.49 |

Every point agrees within |z| ≤ 2.1 on an effect running from 16% down to 2%. The effect is
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
  −1.29% / −0.76% / −0.50% / +0.09% (z −4.1 / −3.9 / −3.1 / +0.7) at 238.6/583/911/2614.5
  keV and FEP within 0.2% (|z| ≤ 0.6). Before source K fluorescence was emitted the totals
  read −22% / −4.1% / −2.6%. The ~1% left at low energy is plausibly Th L x-rays escaping
  the bare sphere's skin (CeeLo emits no L lines; TODO.md).

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

- **HPGe at 59.5 keV behind a thick dead layer is a data difference.** Behind a
  0.0739 cm front dead layer (τ ≈ 0.82) CeeLo reads ~1.0–1.1% below GEANT4 in FEP and
  total at 59.5 keV with either GEANT4 list. GEANT4's photoelectric cross section for light
  and mid-Z elements carries a smooth deficit near 60 keV: CeeLo / GEANT4 at 59.54 keV is
  Ge +1.25%, Al +2.1%, C/O +2.2% (I, Cd, Te, Pb within 0.15%). Against NIST (SRD 126),
  μ/ρ(Ge) is CeeLo +0.22…+0.43% at 30–100 keV (flat) and GEANT4 +0.9% at 30–40 keV,
  −0.72% at 60 keV: CeeLo is the closer. Through the dead layer and endcap that predicts
  −0.97% in transmission.

- **Pair production near threshold.** EPICS2023's pair cross section is 3–9% below
  GEANT4's parameterization at 1332 keV (Ge −9%, I −3%, Na −6%) and 21% below at 1.2 MeV
  for iodine, converging above ~2 MeV (Pb +2.3%, Ge −1.0% at 2614 keV). Single and double
  escape peaks at 1332 keV read 5–13% low against GEANT4 for that reason. Unarbitrated
  (TODO.md).

- **EADL K line energies** come from subshell binding-energy differences: within 0.04 keV
  of the measured lines for Ge, Cd, Te and I, but Pb Kα1 sits at 75.25 keV against
  74.97 keV measured (GEANT4 puts it at 74.97). Visible as a 0.3 keV shift of shield
  fluorescence peaks on HPGe (TODO.md).

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
- **Code placement moves timings by ~10%**: the same library code linked at different
  addresses ran `benchmark_mc_configs` 8–11% slower in CPU time on Apple silicon whenever
  the electron walk landed differently (a driver that only parsed three new flags did
  it). Before attributing a cost change to a code change, rebuild the comparison from the
  same driver and check the hot functions' addresses (`nm`); measure CPU time, interleaved,
  at fixed event counts.
- **Importance-sampling denominator**: `(1/N) * sum(w_i * I_i)`. Using `sum(w_i)` gives
  intrinsic efficiency (it cancels the cone weight), not absolute efficiency.
- **Coincident GEANT4 boolean surfaces** stall navigation ("track stuck") and — worse —
  inflate the `--entry-diag` electron-entry rate ~1.75× while leaving totals unchanged, so
  they silently corrupt channel-level comparisons. The nested source export avoids them;
  keep it that way. (~1.5k benign warnings per run at `active_crystal_PV` remain in every
  config.)
