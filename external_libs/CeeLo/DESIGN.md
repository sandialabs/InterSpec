# CeeLo — design and validation notes

Developer-facing record of how CeeLo is built, what physics it models, how it is
validated against GEANT4, and what is deliberately approximated. For an
introduction to the library see [README.md](README.md); for the list of open
shortcomings see [TODO.md](TODO.md).

## Purpose

CeeLo is a C++17 Monte Carlo photon transport library for computing gamma-ray and X-ray detector efficiency. Given a detector geometry (NaI, LaBr3, CZT, etc.), optional attenuators/shielding, and a source configuration, it computes:

- **Full-energy peak (FEP) efficiency**: probability that a photon deposits 100% of its energy in the detector crystal
- **Total efficiency**: probability that a photon deposits any energy in the crystal
- **Pulse-height spectrum**: energy deposit histogram (compatible with spectral analysis software)
- **Cascade summing corrections**: coincidence summing-in and summing-out factors. Optional
  internal-conversion (IC) electron + coincident K-Auger energy deposition is available via
  `CascadeConfig::enable_ic_electrons` (default OFF, both estimators). It transports the conversion
  electron through the source geometry + air gap and deposits it in the crystal only when it can
  physically reach it (auto distance/solid-angle/air gate), so distant sources are unaffected and the
  off path is bit-identical. Matters at contact/well/4π (e.g. Ba-133 356 keV summing factor drops
  ~0.44→0.31 at contact). GEANT4-validated: flag-ON matches G4 per-decay peak areas to <1% at 2 cm
  (Ba-133), measured with the dev-only IC-electron probe. Two G4 gotchas for the
  cascade cross-check: use an **air world** (`vacuum_world=false`, so G4's IC electrons air-degrade like
  the engine's) and set `/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year` (else a long-lived
  ion parent never decays). Cs-137's Ba-137m isomer does not fully decay under these G4 settings.

The goal is for this library to achieve agreement with GEANT4 11.4.0 reference simulations within ~1% across validated configurations (30 keV to 3 MeV).  Full-energy peak and total efficiencies are the priority; capturing the spectrum accurately is also wanted, and is worth a few percent of slowdown but not much more.


## Build

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j8
cd build && ctest --output-on-failure
```

Tests take ~15-30 minutes (MC simulations with statistical convergence checks).

Dependencies: C++17 compiler, CMake 3.15+, Boost (for Boost.Test), Eigen3.

**Compile-time options** (CMake cache vars): `CEELO_WITH_SANDIADECAY=ON` (cascade adapter),
`CEELO_SOURCE_ESCAPE_MODEL=gate|tails|gs` (source-electron skin-escape model; default `gate` — the
validated empirical gate; `tails`/`gs` are first-principles A/B variants — see Known Limitations →
"Source-electron skin escape" for the per-config accuracy comparison),
`CEELO_BUILD_RESPONSE=ON` (storable detector-response payload — needs header-only **rapidxml**).

**Portability notes** (both learned getting the vendored copy to build with GCC 12 on Linux —
neither shows up on macOS/libc++, so do not "clean them up" without a Linux build):

- The library links `Threads::Threads` **PUBLIC**. `EfficiencyCalculator.cpp` and
  `ResponseGenerator.cpp` use `std::thread`, which on libstdc++ needs `-pthread` at link
  time; on macOS pthread lives in libSystem, so the omission is invisible there.
- RNGs are seeded through the `seeded_rng()` helper in `EfficiencyCalculator.cpp`, which
  uses `std::seed_seq`'s **iterator** constructor. libstdc++ rejects the `initializer_list`
  form for `uint64_t` arguments. Its three words are chosen to reproduce the streams the
  `initializer_list` form produced, bit for bit — adding or reordering words would move
  every seeded result in the library, including the GEANT4-gated validation numbers.

**rapidxml** (only for `CEELO_BUILD_RESPONSE`): NOT vendored — for ODR safety the host
application's own copy is used. The CMake configure searches a few sibling `3rdparty`
directories for `rapidxml/rapidxml.hpp` (SandiaDecay bundles a copy, as does SpecUtils);
when none resolves, pass `-DCEELO_RAPIDXML_DIR=<dir containing rapidxml/>`.

## Open items

[TODO.md](TODO.md) is the single source of truth for known shortcomings, deferred
features, and fidelity/speed/data improvements. Each entry records its **Effect**, what
its **Implementation** would take, and the regime where it **Shows up**; items are
deleted as they are resolved. This document records what *is* modelled and what is
deliberately approximated — not the open-items list.

## Repository layout

The tracked file set is the reusable library: `git archive HEAD` produces exactly the
payload that gets vendored into a host application (library + tests + examples + the
GEANT4 harness + the GEANT4/CeeLo reference results). CeeLo's own efficiency results
live in `tests/data/ceelo_reference/`, paired 1:1 with the GEANT4 references in
`tests/data/geant4_reference/`; `profiling/compare_validation.py` gates them against
each other.

One-off studies, parameter sweeps, external-code benchmarks, and design notes live in a
`studies/` tree that is **not** part of the distributed source. The root `CMakeLists.txt`
adds it only `if(EXISTS studies/CMakeLists.txt)`, so it is a no-op here. A few comments
in the source refer to those studies by name for provenance; the conclusions they reached
are summarised in this document.

The cascade adapter's SandiaDecay dependency is located via the `CEELO_SANDIADECAY_DIR`
CMake variable — point it at the host project's copy when vendoring.

## Source Layout

```
src/
  cascade/          True-coincidence (cascade) summing: analytic + correlated-MC
                    estimators, beta/IC spectra, optional SandiaDecay adapter
  cross_sections/   Compiled-in photon cross-section data (auto-generated
                    element_data.cpp) and its lookup/interpolation
  efficiency/       EfficiencyCalculator - main entry point
  export/           GDML and macro export for GEANT4 validation
  geometry/         Ray tracing (cylinder, box, sphere), attenuator shells, source geometry
  io/               Storable detector-response payload, MC-free response kernel,
                    EFFTRAN-style efficiency transfer (CEELO_BUILD_RESPONSE)
  materials/        Material definitions (composition, density), builtin materials
  physics/          Electron CSDA range tables, stopping powers, Moliere scattering
  transport/        Photon transport (Compton, PE, Rayleigh, PP, fluorescence, bremsstrahlung)
tests/              Boost.Test suites (~35 ctest targets)
tools/
  prepare_cross_sections/  Locked EPQ/EPICS2023/ESTAR generators and reports
  geant4_validation/    G4 harness and aggregate-observable comparison tools
visualization/      Dependency-free WebGL viewer for exported GDML geometries
```

## Key Design Decisions

### Physics

- **Analog MC transport**: No variance reduction inside the detector crystal. Each photon interaction (Compton, PE, Rayleigh, pair production) is sampled explicitly. Electron energy deposition uses CSDA (continuous slowing-down approximation) with Moliere multiple scattering for KE > 200 keV.

- **Cone importance sampling**: For point sources with no shielding, photons are emitted within a cone subtending the detector (not isotropic 4pi). Each event carries a weight = cone solid angle fraction. This gives ~100x speedup for typical geometries. The standard IS estimator is `(1/N) * sum(w_i * I_i)` -- the denominator is N (number of events), NOT sum of weights.

- **Cross-section data**: Photon magnitudes and angular factors are generated
  directly from locked **EPICS2023 EPDL** records. Each process has its own
  per-element grid with source nodes, thresholds, and both sides of absorption
  discontinuities retained. Adaptive interval probes gate maximum material
  representation error at 0.1% and integrated error at 0.05%; angular tables
  additionally gate CDF and mean-cosine errors. Generated tables must not be
  hand-edited; regenerate and audit them through `tools/prepare_cross_sections/`.

- **K-shell fluorescence**: After photoelectric absorption, K-shell vs
  outer-shell is sampled using direct EPICS2023 EPDL MT=534/MT=522 fractions.
  K-shell vacancies produce fluorescence X-rays or Auger electrons using direct
  EPICS2023 EADL MF=28/MT=533 relaxation data through Z=99, covering
  radioactive-decay daughters beyond the Z=92 transport-table limit.

- **Bremsstrahlung**: Seltzer–Berger scaled differential cross-sections (92
  elements, **27 native electron energies from 10 keV–20 MeV × 32 κ**).
  Rejection sampling uses a 1/k envelope and transports brems photons
  recursively. The `float32` production table comes from pinned public-domain
  NIST EPQ. A fixed transformation constructs κ=0.025 quadratically through
  EPQ {0, 0.05, 0.1} and κ=0.075 cubically through {0, 0.05, 0.1, 0.15}; all
  remaining values are copied from equal-κ EPQ nodes. There is no row
  normalization or element/energy-specific correction. This aligns the low-κ
  grid with the current Geant4 grid without extracting or fitting G4EMLOW
  coefficients. Aggregate Geant4 spectra/rates are retained only as validation
  observables, with manual shape approval required before feature migration.

- **Compton scattering**: Butcher-Messel method (same as GEANT4 G4LowEPComptonModel) for the scattering angle, with S(x,Z) bound-electron rejection. A previous implementation using Kahn's method had a subtle bug that caused +5-7% FEP overestimates. The scattered-photon energy is Doppler-broadened by the bound electron's momentum (impulse approximation, analytic subshell profiles) — see Known Limitations → Doppler broadening.

### Distances, frames and reference points

Every distance in CeeLo is one of three things, and confusing them is the single
easiest way to introduce a silent few-percent bias. They do not crash and no
test of internal self-consistency catches them — only a comparison against an
external truth does.

- **The crystal frame.** The one frame the tracer works in. `z = 0` is the front
  face of the **crystal solid**, `z = L` its back face; a source in front is at
  **negative z**. Radial distance is from the symmetry axis.

- **The dead layer is carved out of the inside of the crystal.** With a dead
  layer the solid still spans `[0, L]` at radius `R`; only the *active* volume
  shrinks, to `[t_front, L - t_back]` and radius `R - t_side`. Two consequences
  that have each already caused a bug: the crystal face stays at `z = 0` no
  matter how thick the dead layer is, and **a dead layer must never be added to
  an outward distance**. "The front of the crystal" and "the front of the active
  volume" are different planes `t_front` apart — always say which you mean.

- **Only attenuator shells sit in front of the crystal.** Each extends its front
  thickness toward negative z, so the outermost front surface — the "endcap
  front", the face a user can touch — is at `z = -Σ(layer front thicknesses)`.

`ReferencePoint` selects which plane a caller's `dist_cm` is measured from, and
`GeometryDescriptor::endcap_front_offset_cm()` is the **only** conversion
between them, applied in `DetectorResponse::query_position()`. Anything else
building a source position from a user distance must go through that function or
reproduce it exactly.

That function used to add the dead layer to the offset, which placed every
`EndcapFront`-referenced source one dead-layer thickness too far from the
crystal — ~0.7 mm on a typical HPGe. In the far field that is invisible; at
contact geometry it is a near-uniform ~3 % efficiency deficit above ~100 keV,
which is exactly what it looked like when finally caught against ANGLE's own
numbers for a GEM35-70 (`EFFTRAN` vs ANGLE went from a 2.5 % median error to
0.9 % when it was fixed). Below ~50 keV it does not show up at all: the crystal
is effectively large enough there that its solid angle has saturated and is
insensitive to a fraction of a millimetre.

Two nearby functions look like the same bug and are **not** — do not "fix" them:

| function | includes dead layer? | why |
|---|---|---|
| `endcap_front_offset_cm()` | **no** (fixed) | a real distance, query-time only, so correcting it costs nothing |
| `transverse_half_extent()` | yes, deliberately | a *scale parameter*, not a dimension; computed identically at generation and query time, so changing it desynchronises stored responses |
| `Geometry::outer_bounding_radius()` | yes, deliberately | a sampling *bound*; over-estimating is safe, under-estimating loses paths |

The rule of thumb: if a quantity is used at generation **and** query it is a
model parameter and must not change; if it is used only at query time it is an
interpretation and can be corrected.

### Geometry

- **Detector shapes**: Cylinder and box. Defined by dimensions (radius+length or x/y/z half-widths). Crystal front face at z=0, extends along +z.

- **Bulletized (rounded) front edge**: HPGe crystals are usually *bulletized* -- the outer front edge is a quarter-torus fillet rather than a sharp 90° corner (ANGLE `.outx` files call it `bulletizingRadius`). `set_bullet_radius(r_b)` turns it on for cylinders; 0 (the default) is a sharp edge. The fillet's ring circle has radius `rho_c = R - r_b` and lies in the plane `z = r_b`, tangent to the front face and to the side wall. Removed material, by Pappus: `V = 2 pi rho_c r_b^2 (1 - pi/4) + pi r_b^3 / 3` (1.28% of the crystal for the GEM35-70 example: R=2.915, L=6.89, r_b=0.8 cm).

  The solid is **convex and a subset of the sharp cylinder**, and `Cylinder.cpp` exploits that rather than intersecting a torus outright: the sharp interval is computed exactly as before, and only an endpoint that lands *inside the removed corner* is moved onto the fillet. Endpoints on a shared surface (front face inside `rho_c`, side wall above `z_c`, back face) are already exact and are left untouched bit for bit. The fillet root is found by safeguarded Newton on `G(t) = (rho - rho_c)^2 + (z - z_c)^2 - r_b^2`, which is **convex wherever rho >= rho_c** -- i.e. on every corner-zone sub-segment -- and whose bracket is guaranteed because every corner-zone exit point inside the sharp cylinder is within `r_b` of the ring centre. Chords that stay in the corner zone end to end (grazing corner clips) are resolved by minimising the same convex `G`. The closed-form ray-torus quartic is deliberately **not** used: its coefficients cancel catastrophically for near-double roots, which is exactly the grazing-corner case. It is kept as a comment and as an independent test oracle.

  Because every consumer -- photon transport, the electron CSDA/Molière walk (which re-traces per substep), `ResponseKernel`, the stochastic-Rayleigh re-trace -- goes through `Geometry::trace_ray`, implementing the shape in `trace_cylinder_geometry` covers all of them. Cone biasing and solid angle are untouched: they use the bounding cylinder, which bulletization only shrinks material inside of. Attenuator shells stay sharp, which is physically right -- the endcap does not follow the crystal fillet.

  With a dead layer, the active volume gets the inward-offset fillet: radius `r_b - max(t_front, t_side)` about a centre placed to stay tangent to both offset faces. For equal front/side thicknesses (the usual case) that is the exact offset surface -- same ring centre, radius `r_b - t`. Dead layers are ~0.1 mm against an ~8 mm fillet, so the convention for unequal thicknesses is far below any measurable effect.

  **Performance**: sharp crystals pay one predicate (`bullet_radius_ > 0.0`) and take the original code paths, verified bit-identical on fixed-seed runs of a bare-NaI and a bored-HPGe geometry. Timing on config 1 (2M events, single thread, min of 3) went 9.03 s -> 8.20 s across the change: no regression -- the small speedup is a code-layout side effect, since the results are unchanged bit for bit.

- **Rounded bore tip**: `set_bore_hole(r, depth, rounded_tip=true)` caps the bore's closed end with a hemisphere of the bore radius (a round-tipped drill; ANGLE's `<core rounded="yes"/>`). The stated depth is preserved -- the hemisphere's apex sits where the flat bottom would -- so rounding returns `pi r^3 / 3` of crystal (0.13 cm³, ~0.07%, for the GEM35-70 bore). The capped bore is convex and its upper half lies inside the straight part, so the ray interval is just the hull of the straight tube's and the ball's intervals: a closed-form quadratic, no iteration.

- **Stored geometry (`GeometryDescriptor`)**: the serializable form of all of the above, and the only thing a saved response carries -- `build_geometry()` (and `ResponseGenerator::configure_calculator()`) replay it onto a `Geometry` / `EfficiencyCalculator`. It holds `bullet_radius_cm` alongside `dimensions_cm`, and the rounded-tip flag inside its `BoreHoleConfig`. Both are written to XML only when non-default (`bulletRadius` on `<Detector>`, `roundedTip` on `<Bore>`), so responses for sharp-edged, flat-bored crystals keep their exact bytes and `content_hash`; a reader that predates the attributes, or a file that predates them, defaults to sharp-and-flat, which is what those crystals actually were. `sm_xmlSerializationVersion` is deliberately **not** bumped for such purely additive attributes -- the only version logic is a range check, so a higher version would be rejected outright by existing builds instead of degrading gracefully.

  `Geometry`'s preconditions (fillet vs radius/length/dead layer, bore vs radius/depth/fillet) are `assert`s, so a bad descriptor would trace silent garbage in a release build. `GeometryDescriptor::problems()` is the release-safe mirror of those asserts, plus the cases they miss: a `dimensions_cm` too short for the shape (`set_detector` only asserts its length, then indexes it), a dead layer that consumes the crystal, a bore wider than the *active* radius `R - t_side`, and a bore that clears the outer fillet but not the dead-layer-offset *active* fillet. The last two are policy as much as mirroring -- they produce zero active volume rather than tripping an assert, and refusing to load beats silently computing zero efficiency. `build_geometry()` and `configure_calculator()` both run it and throw. Callers building geometry from user input or an imported file should consult it first and report or relax rather than throw (InterSpec's ANGLE import drops an unrepresentable fillet or rounded tip, with a warning, instead of failing the import).

- **Attenuators** (detector shielding): Cup-shaped shells around the crystal, added via `add_attenuator(material, front_thickness, side_thickness, z_start, z_end)`. These are concentric -- each layer wraps around the previous. Attenuators are part of the detector geometry and are always traversed in ray tracing.

- **Source shielding** (around the source): Spherical shells (point source) or concentric shells matching source shape. Added via `add_source_shield(material, thickness)` (uniform), `add_source_shield(material, t_radial, t_end)` (cylindrical sources), or `add_source_shield(material, t_x, t_y, t_z)` (rectangular sources). Per-dimension thicknesses apply symmetrically to both ± faces of each axis; individual components may be zero (but not all), e.g. `(t_radial, 0)` = side wall with open ends. Call `set_*_source()` before adding per-dimension shields. Source shielding uses full MC transport in full mode (Compton/PE/Rayleigh with re-tracing after scatter) and stochastic Rayleigh in FEP-only mode.

  Full-mode source transport also generates **secondaries** for all source shapes (`SourceFullTransportResult::secondaries`, processed by the detector-transport loop in `simulate_thread`): (1) 511 keV annihilation gammas from pair production in source material/shields; (2) bremsstrahlung from Compton-recoil / photo / pair electrons stopping in the source geometry (Seltzer-Berger sampling along a Molière condensed-history walk, `ElectronCsda::walk_in_source_geometry`); and (3) escaped electrons that cross the gap into the crystal (Molière walk escape detection; enabled per-config via `enable_source_electron_transport`). The Molière direction diffusion is essential for the electron-escape channel — straight-line CSDA overestimates escape from dense shields (Fe) by 2.4–10x. **All source shapes (incl. Marinelli) now use this walk** for electron escape — the legacy straight-line `try_source_electron` was retired when Marinelli migrated to `walk_in_source_geometry` (water → the regime-aware exit-state gate runs fully analog), which let config 8 join the gated set; `min_distance_to_boundary` has a Marinelli case so the containment fast-path applies there too. These channels matter only above the 1022 keV PP threshold and/or for recoil electrons > 200 keV; for config 11 (0.5 cm Fe) at 3 MeV they contribute ~+4.2% (annihilation), ~+3% (brems), and ~+1% (electron escape) of total efficiency. A `set_source_brems(false)` toggle and `benchmark_mc_configs --no-source-brems` exist for A/B quantification. **Perf**: the geometry-tracing Molière walk (it re-traces the source segments every substep) dominated these channels' cost. `process_electron` now takes a cheap path for electrons that stop inside their birth material — the common case: an exact containment test (the electron's CSDA range in its birth material `mat`, `range/density`, < distance to the `mat`-region surface, `SourceGeometry::min_distance_to_boundary`) routes them to the geometry-free `sample_brems_in_material` (same SB physics, creation-point emission — the Marinelli sampler) and skips escape (provably impossible). Only electrons that can reach an interface or escape take the full walk. ~1.2–1.25x faster on configs 11/12 above 1 MeV; totals unchanged within statistics. Further speedups would have to touch the genuine transport of escape-capable electrons and emitted brems photons.

- **GDML export**: Attenuators are exported as cup-shaped subtraction solids (outer full cylinder minus shorter inner cylinder) to properly cover the crystal front face. This avoids the hollow-tube bug that caused -26% discrepancy at 30 keV.

  A crystal with a **bore and/or a bulletized front edge** is a surface of revolution, so it is exported as one native `<polycone>` named `CrystalSolid` rather than a stack of boolean solids: the `rmax` profile follows the fillet arc up to `R` over the first `r_b` of depth, and the `rmin` profile is 0 until the bore starts and then the bore radius (via the tip hemisphere first, if round-tipped). Circular arcs become 64 chord segments, which keeps the sagitta well under a micron (~1e-4 % of the crystal volume — verified against the tracer's own volume across bore/fillet combinations). Plain cylinders still export as a single `<tube>`, so every pre-existing config's GDML is byte-for-byte what it was.

  This is not just tidiness. An earlier subtraction-based version of this export (tube minus a corner-ring-minus-torus cutter, plus a tube∪orb bore) produced **~3250 stuck-track warnings per million** under cone-biased runs — GEANT4 navigation struggles with tracks skimming a subtracted surface — which both floods the log and risks biasing the very comparison the export exists to support. The polycone produces ~25 per million, the same rate as a plain tube. Two consequences worth remembering: the two profiles must be merged onto the *union* of their z values rather than truncating one at the other's start (otherwise a bore whose closed end reaches into the front fillet silently cuts the arc short), and the harness's `--cone-bias` bounds parser has to understand `<polycone>` as well as `<tube>`.

  **The bore cavity must not be filled with a daughter volume.** The bore is already subtracted from the crystal solid, so an `Air` daughter placed in it lies *outside* its mother and G4 navigation gets stuck on it. The cavity is left to the world material, which matches the MC (the bore region gets no material there either) whenever the export uses `vacuum_world=true`, as the validation configs do. This bug was latent until configs 25/26: no earlier exported config used `set_bore_hole`. (Note also that GDML resolves `volumeref` only against volumes already read, so any daughter must be *defined before* the volume that places it.)

### Simulation Modes

- **Full mode** (default): Tracks all photon interactions through attenuators and crystal. Electron CSDA with bremsstrahlung. Gives both FEP and total efficiency, plus pulse-height spectrum.

- **FEP-only mode**: Full MC transport inside the scoring crystal (same Compton/PE/fluorescence physics as full mode), but uses importance-weighted stochastic Rayleigh through non-scoring segments (attenuators, source shielding) -- treating Compton/PE as absorption and only sampling Rayleigh scattering. Photons that exit the crystal boundary are immediately killed (can't contribute to FEP). Faster than full mode but only computes FEP efficiency (no spectrum or total efficiency).

### Efficiency transfer (EFFTRAN-style, `CEELO_BUILD_RESPONSE`)

`src/io/EfficiencyTransfer.h` packages the classic EFFTRAN workflow: anchor an
efficiency curve at ONE reference position (from CeeLo MC **or** user-measured
`(E, eff)` points) and transfer it to any other position via the ratio of the
interaction-weighted effective solid angle,
`eps(E, target) = eps_ref(E) · K(E, target) / K(E, ref)`, where
`K = ApertureQuadrature::interaction_omega` (the same MC-free
`ResponseKernel` the parameterized DRF already uses). This is EFFTRAN's `T`-ratio
transfer, but fully 3-D (off-axis native), Rayleigh-aware, and chord-exact —
holding the response factorization's residual `eta(E) = eps_ref/K` constant
across positions. Two entry points: the standalone `EfficiencyTransfer` class
(live cross sections; `eps_at`, `eps_at_with_tsrc` for MEFFTRAN extended-source
self-absorption), and `make_transfer_response(...)`, which emits a normal
`DetectorResponse` (angle-flat eta over stored mu tables) so XML / extended-source
/ cascade / grounding all work unchanged. `ResponseGenerator` gains a
`GenerationOptions::transfer_mode` (+ `n_anchor_angles`) that builds a DRF from
MC at only the on-axis backbone (~8 nodes flat; a few forced cos-theta anchors
when `n_anchor_angles > 1`), vs. hundreds for the General profile. **Validity**:
the K ratio transfers distance at fixed angle to ~1% for d ≳ a few `a` and to
<0.3% far-field (campaign D1); it has no near-field N model (below ~2a it flags
`NearFieldUnmodeled`) and a single on-axis anchor does not carry the `eta(E,theta)`
angular residual (campaign D2). Both regimes are covered by an unconditionally-applied
`DetectorResponse::model_transfer` (`SigmaTransferModel`) that inflates the
reported sigma off-axis/near. `tools/efftran_transfer_demo` reproduces EFFTRAN's
`et.bat`: MC once, transfer over a distance×angle ladder, compare to direct MC
with z-scores. Tests: `tests/test_efficiency_transfer.cpp` (MC-free identities +
a `CEELO_RUN_MC_TESTS`-gated consistency case).

**Measured accuracy/timing** (July 2026; harness + report in the dev-only
the dev-only transfer study; NaI 3×3 and a representative HPGe coax — 0.5 cm bore, 0.7 mm
Li dead layer, 1.5 mm Al endcap + 5 mm gap — transferring a far-field on-axis
anchor vs direct MC at ~0.6% stats): **FEP** on-axis ≤1% for d≥5a (median
0.2–0.55%, p95 ≤1.8%), ~1.5–2% at 3–5a, ~3–3.6% median (≤5%) near contact (d<3a);
**off-axis** ~3.2–3.6% median / ~6% p95 from a single on-axis anchor, dropping to
**<1% median** (p95 2–3%) with 3 forced cos-theta anchors (`n_anchor_angles=3`,
~+25 s MC). **Total efficiency**: on-axis ≤0.6%, NaI ~1% everywhere. The HPGe
*near-field* total under-prediction (dead-layer/endcap Compton scatter-in the
removal-only kernel misses; −2.5% at contact) is corrected by a **passive-layer
scatter-in credit** in the total kernel (`interaction_omega`'s
`passive_compton_recapture`; `DetectorResponse::scatter_in_recapture`, default
`kTotalScatterInRecapture = 0.8`, calibrated in the transfer study): it folds
`(mu_total − r·mu_cs)` on non-scoring segments so the credit grows with the
oblique/near-field passive path and vanishes far-field — HPGe contact total
−2.5%→−1.4%, NaI unaffected, FEP bit-identical (peak photons can't be
scatter-in). HPGe total *off-axis* still drifts (the angle-flat tier can't carry
the angular η_tot shape) — use the few-anchor mode's η_tot table there. **Timing** (matched 0.01 per-node precision, auto-threaded): a
transfer-flat DRF is ~45 s / 28 MC nodes vs a full General-profile DRF at 496 s
(NaI) / 714 s (HPGe) / 544 nodes — **11–16× faster**; a user-supplied *measured*
curve needs **zero** MC (sub-second deterministic transfer). Bottom line for a
host (InterSpec): viable as a "quick DRF" or a measured-curve → other-geometry
estimator; on-axis-far is ≤1%, near/off-axis is a few % with honest inflated σ +
provenance flags, and the full MC stays available as an upgrade.

## Cone Biasing Policy

**Use pure cone biasing only for point sources with no source shielding or source material — or when using FEP-only mode.**

When source shielding or extended sources (non-zero density) are present in full mode, photons can scatter in the source material and change direction. Pure cone biasing would miss these scattered photons since they may originate from directions outside the cone. In full mode with source effects, emit isotropically into 4π — or use **mixture angular biasing** (`BiasingConfig::mixture_cone_alpha`, see Event Biasing below), which is unbiased even with source scatter because its emission pdf q(Ω) = α·(1/4π) + (1−α)·(1/Ω_cone) is nonzero over all 4π.

In FEP-only mode, cone biasing is safe even with source material, because the stochastic Rayleigh transmission weighting handles source effects without changing direction (most of the time). Per-position cone sizing is used for extended sources: each event computes its own cone from the sampled source position, giving much tighter cones for positions near the detector axis.

This applies to both our MC code and the GEANT4 validation harness:
- `generate_all_spectra.cpp`: `run_config(..., use_cone_bias)` controls both MC histogram generation and the generated `run_g4_all.sh` script
- G4 harness: `--cone-bias` flag on the command line
- Isotropic G4 runs need ~50-100x more events to achieve comparable statistics in the detector

## Event Biasing (Variance Reduction)

`BiasingConfig` (set via `EfficiencyCalculator::set_biasing()`) provides additional unbiased
variance-reduction techniques beyond cone/position biasing. When `set_biasing()` is NOT called,
`compute()` auto-enables them per geometry, energy, and termination targets
(`compute_effective_biasing()`); pass a default-constructed `BiasingConfig` to force fully
analog sampling (used by A/B consistency tests).

- **Forced first interaction** (`force_detector_interaction`): the primary photon's first
  interaction is forced inside the detector geometry; the history carries weight
  `1 − exp(−τ)` ≤ 1 (`TransportResult::weight`), with the interaction point drawn from the
  truncated exponential across path segments. Secondaries are never re-forced (the flag is
  cleared in transport's local config copy before any recursion). Gated off for the Marinelli
  primary path (re-entry physics needs the un-collided pass-through branch). For source-effect
  events the gate is **per event** (`force_this_event` in `simulate_thread`): forcing is valid
  only when the primary is the sole depositing branch, so events whose source transport
  produced secondaries (annihilation gammas, shield brems) or escaped electrons run the
  primary analog. Unbiased because the estimator choice conditions only on already-sampled
  source randomness.
  **Measured**: wins where pass-through probability is high — CZT ≥ 200 keV (FEP FOM ×1.2–1.4,
  total ×1.5–2.7); *hurts* FEP FOM ×0.75–0.95 in thick crystals (every event pays full
  transport cost) while helping total ×1.3–1.8. Auto-enable: central-ray optical depth
  τ < 0.7 when FEP precision matters; τ < 4 when only total efficiency is targeted.

- **Mixture angular biasing** (`mixture_cone_alpha` = α): for full-mode runs with source
  effects (non-Marinelli), emit from q(Ω) = α·(1/4π) + (1−α)·(1/Ω_cone); weight
  w = (1/4π)/q(Ω) ≤ 1/α. Replaces pure-isotropic emission; unbiased for source-scattered
  photons. Isotropic-branch draws landing inside the cone must get the in-cone weight (the
  membership test shares the cone sampler's axis, `normalize(−src_pos)`). Degenerate cones
  (half-angle ≥ π − 0.01) fall back to isotropic with w = 1.
  **Measured**: FEP FOM ×5.6–6.8 (config 11) and ×17–29 (config 12) at α = 0.1–0.2, but total
  efficiency FOM *degrades* (×0.04–0.6) because the scattered-in component is sampled only
  α of the time at weight 1/α. **No longer auto-enabled** — superseded by the two-stream
  estimator (below); kept as a manual `set_biasing()` knob.

- **Two-stream stratified estimator** (`two_stream`, `two_stream_direct_fraction` = f): for
  full-mode source-effect runs (non-Marinelli), every tally partitions exactly by the
  primary's source-transport outcome (`SourceFullTransportResult::interacted`): the u-term
  (zero interactions — straight line, full energy, no secondaries; the only meaningful FEP
  contributor plus the direct part of total) and the s-term (interacted — the wide-angle
  scattered-in component that dominates total under shielding; it still carries real FEP:
  small-angle Compton at ≤100 keV and Rayleigh-only chains stay inside the ±1.5 keV window,
  up to ~24% of FEP for cellulose at 59 keV — never assume FEP_s = 0). A deterministic
  round-robin fraction f of events runs the **direct stream**: cone emission weighted by
  (Ω/4π)·T, T = exp(−Σ μ_tot·ℓ) = `SourceGeometry::no_interaction_probability()` (exactly
  the analog zero-interaction probability — source transport samples unmasked μ_total), then
  plain detector transport at E₀ (forcing always valid: provably sole branch). The rest run
  the **scatter stream**: isotropic analog source transport, killed early when the primary
  did not interact (required for unbiasedness of the split + a time win), else the unchanged
  legacy detector processing. Stream weight factors 1/f and 1/(1−f) keep the single-pool
  estimator Σw·I/N unbiased for FEP/total/spectrum with no tally-machinery changes; the
  reported variance is slightly conservative (between-strata term of the iid formula).
  **Measured** (vs same-session analog): at f=0.25 + Compton bias γ=0.3, joint
  min(FOM_fep, FOM_tot) ×2.9–6.1 (cfg 11), ×5.4–15.7 (cfg 12), ×5.9–17.3 (2 cm Pb stress),
  with total FOM ×1.8–3.8 (the mixture's total regression is reversed); FEP-only-target
  f=0.5 keeps FEP FOM ×5–30 (×670+ for tiny-detector/big-source). f ≥ 0.85 is dangerous at
  low energy (FEP_s undersampling: cfg-12-like @59 keV FEP FOM ×0.16). Auto-enable: f=0.5
  when only FEP precision targeted, else f=0.25 — but only when the detector cone from the
  source center subtends Ω/4π < 0.15 (close geometries like a puck 1 cm from the face
  measurably lose and stay analog).

- **Compton-angle mixture biasing** (`compton_cone_fraction` = γ): at the primary's FIRST
  Compton vertex in the source geometry, the outgoing direction is sampled from
  q(ω) = (1−γ)·p(ω) + γ·1[ω ∈ cone]/Ω_d, where p = KN×S(x,Z) density on the sphere and the
  cone subtends the detector bounding sphere from the vertex; the event weight gains
  p/q ≤ 1/(1−γ). Energy/recoil at the sampled angle use the exact analog kinematics
  (`finish_compton_at_angle`; analog RNG draw order unchanged). The angular norm is composite
  GL32 in s = sin(θ/2) (the substitution removes the √(1−μ) cusp S(x,Z) imprints at μ→1;
  ~1e-6 relative, unit-tested) and cached — first-vertex-only keeps every norm call at the
  constant primary energy. Annihilation/brems recursion and Marinelli are never biased.
  **Measured**: total FOM ×1.8–3.8 on shielded configs on top of two-stream; γ ≥ 0.5
  degrades large unshielded sources (20 cm water: total ×0.52 at γ=0.7). Auto-enable: γ=0.3
  iff total/spectrum requested (γ only helps totals).

Weight-composition rules (apply to any new technique):
- One weight per event: every recursive `transport_photon()` call belongs to the same history;
  the weight is applied exactly once at the tally.
- The estimator denominator stays N (`sum_weights += 1.0` per event) — see Common Pitfalls.
- Branch-level Russian roulette is **fundamentally unsafe for FEP**: the FEP indicator
  I(Σ deposits = E₀) is non-linear in per-branch deposits, so killing one branch with survivor
  reweighting biases FEP even though it preserves the mean total deposit. Only whole-event
  roulette/splitting is valid.
- FOM benchmark harness: `build/examples/benchmark_mc_configs --fom-csv <file>` with
  `--analog` / `--force-interaction` / `--mixture-alpha <a>` / `--two-stream <f>` /
  `--compton-bias <g>` / `--precision <p>` for A/B comparisons (FOM = 1/(σ_rel² · T));
  physics-channel A/B via `--no-source-brems` / `--no-source-electrons`; precision-targeted
  runs also write `our_<cfg>_multi.csv` (consumed by `profiling/compare_validation.py`);
  consistency tests in `tests/test_biasing_consistency.cpp`. Out-of-sample stress configs
  13–17, 19 in `examples/bench_configs.h` (thin Al / 2 cm Pb / 20 cm water / close puck /
  CZT+soil cube / off-axis Fe) guard the auto-policy outside the G4-anchored set.
- Per-event time ceiling (profiled June 2026): the source-electron brems machinery
  (`sample_brems_in_material`, `residual_energy_keV`, `brems_pemit_corr`, Molière walk) is
  40–60% of CPU on source-electron configs; further FOM there is perf work, not estimator work.

## Source Geometry vs Detector Geometry

The simulation distinguishes two separate geometry systems:

- **Detector geometry** (`Geometry` / `RayTrace`): The crystal plus attenuator shells. Ray tracing returns `PathSegment` vectors with material, thickness, and `is_scoring` flag. Both full and FEP-only modes traverse these segments. Attenuators are added via `add_attenuator()` on `EfficiencyCalculator`.

- **Source geometry** (`SourceGeometry`): Source material (the radioactive volume itself) plus source shielding shells around it. Configured via `set_source_material()` and `add_source_shield()`. These are geometrically separate from the detector -- there is open space between the outermost source shield and the detector front face.

### Source shapes (`SourceType` / `SourceGeometry::Shape`)

`Point`, `Cylindrical`, `Rectangular`, `Spherical`, `Marinelli`. All extended shapes self-attenuate (`set_source_material`), take concentric shape-matched shields, and support uniform or exponential-depth activity. Shape-dispatch is centralized: `sample_source_position`, `source_material_path`, `trace_source_segments`, `compute_transmission[_fep_only]`, `min_distance_to_boundary`, GDML export. **Both `compute` and `compute_cascade` pick up every shape for free** (cascade routes through `sample_source_position` + `transport_source_photon`).

- **Spherical** (`set_spherical_source(center, radius, rotation=I, inner_radius=0)`): solid ball or, with `inner_radius > 0`, a hollow **shell** (non-attenuating void center — round flasks, spherical phantoms, void-center shells). Uniform-in-shell sampler `r = cbrt(r_in³+(r_out³−r_in³)U)`; spherical shields reuse the point-source shell tracer (`spherical_shell_path`). Rotation is stored but physically irrelevant (rotation-invariant).
- **Hollow / annular cylinder** (`set_cylindrical_source(..., inner_radius>0)`): tube/pipe/ring with a non-attenuating central bore. Annular sampler `r = √(r_in²+(r_out²−r_in²)U)`; `source_material_path` charges only the annulus (subtracts the bore). `inner_radius = 0` reduces bit-identically to the solid cylinder.
- **Hollow rectangular box shell** (`set_rectangular_source(..., inner_half_dims)`): crate/container walls — outer box minus a centered inner void box (same rotation, non-attenuating). Rejection sampler over the outer-box proposal (acceptance `1 − Vi/Vo`); `source_material_path` subtracts the void chord (`box_shell_path`); GDML fill exports as a box subtraction solid. All-zero `inner_half_dims` reduces bit-identically to the solid box. See `tests/test_hollow_box_source.cpp`; pipe-modeling conventions (contents vs wall, side-on shield composition) in `tests/test_pipe_source.cpp`.
- **CylinderEndOn vs CylinderSideOn** (InterSpec orientations) are set by the `set_cylindrical_source` rotation: **EndOn** = identity (cylinder axis ∥ detector axis, the default); **SideOn** = 90° about y, `R = [[0,0,1],[0,1,0],[-1,0,0]]` (local-z ← detector-x, detector views the curved side).
- **Trace vs self-attenuating** are not separate shapes — see `set_source_material` doc: efficiency depends only on the spatial distribution (uniform / exponential-depth) + self-attenuating medium; InterSpec's per-cm³/g/m² concentration and mass-fraction are downstream activity normalization (Bq×volume) that does not affect efficiency.

The photon's journey: source emission point -> source material -> source shields -> air gap -> detector attenuators -> crystal. In full mode, source geometry uses analog MC transport (`transport_source_photon()`); in FEP-only mode, it uses stochastic Rayleigh weighting (`compute_transmission_fep_only()`). The detector geometry always uses full ray tracing and segment-by-segment transport.

When source geometry has effects (`has_source_effects() == true`) and the simulation is in full mode, photon emission is isotropic 4pi (not cone-biased), since scattering in source material can redirect photons toward the detector from any initial direction.

**Note**: Air between the outermost source shielding layer (or the source itself if no shielding is present) and the detector is **not simulated** in the MC transport. Air attenuation/scatter effects will eventually be applied as a deterministic (analytical) correction after the simulation.

## GEANT4 Validation

### Harness

**What you need:** a GEANT4 11.4.0 installation (`bin/geant4.sh` to set the environment
and data paths before any run; `lib/cmake/Geant4` for the CMake config). It is
required only for validation; production coefficients come from the locked EPQ,
EPICS2023, ESTAR, and xraylib sources. Having the GEANT4 *source* tree available
is useful for comparing model choices (see "Inspecting G4 internals" below).

The harness source is `tools/geant4_validation/` and builds separately against GEANT4
(done once):
```bash
source <your-geant4-install>/bin/geant4.sh
cmake -B <g4val-build> -S tools/geant4_validation \
      -DGeant4_DIR=<your-geant4-install>/lib/cmake/Geant4
cmake --build <g4val-build>
```
(On macOS, if a run aborts with exit 134 / a dyld error: the GEANT4 libraries dynamically
link Qt, so `brew install qtbase`.)

**Running a G4 comparison — full workflow:**
```bash
# 1. Generate the G4 inputs (GDML geometry + GPS macro) from the engine, so the G4 geometry
#    and source EXACTLY match the MC. Exportable configs: 1-7, 12, 20-24 (benchmark_mc_configs);
#    --events tiny just to trigger the export, --macro-events sets the G4 /run/beamOn count.
cd build/examples
./benchmark_mc_configs --config 20 --energies 662,2614 --events 20000 --macro-events 12000000
#    -> writes detector_20.gdml + run_20_662keV.mac + run_20_2614keV.mac
#    (Marinelli/cfg 8 + cfg 11 are exported by generate_all_spectra instead — see "Why some
#     cfgs aren't exportable": can_export in benchmark_mc_configs.cpp is a vetted opt-in list,
#     because an extended source's GPS sampling volume must EXACTLY equal the engine source
#     region or escape/self-absorption is silently biased.)

# 2. Run G4 (MUST source the env first). Point sources: add --cone-bias (matches the MC cone).
#    Extended/shielded sources: ISOTROPIC (no --cone-bias) and need 16M+ events for stats.
source ~/install/GEANT4-11.4.0/bin/geant4.sh
~/install/src/g4val-build/ceelo_g4val detector_20.gdml run_20_662keV.mac g4val_cfg20_662keV.csv

# 3. Compare: assemble a tests/data/geant4_reference/<name>_multi.csv from the per-energy CSVs,
#    add the config to profiling/compare_validation.py CONFIGS, then run it against our_N_multi.csv
#    (written by `benchmark_mc_configs --config N --precision p`).
python3 ../../profiling/compare_validation.py
```
Output CSV columns: `energy_keV,fep_efficiency,fep_uncertainty,total_efficiency,total_uncertainty,num_events`.
Macros carry a `/run/beamOn <N>` line set by `--macro-events`; `sed` it to change the count without
re-exporting.

**Radioactive-ion (cascade) source macros — two gotchas that silently give wrong/empty output:**
For a full-decay cascade source (`/gps/particle ion` + `/gps/ion <Z> <A>`, used for true-coincidence
summing validation, e.g. `56 133` = Ba-133):
- **`/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year`** must be in the macro (right after
  `/run/initialize`) — otherwise a long-lived parent (Ba-133 10.5 yr, Cs-137 30 yr, most calibration
  nuclides) is treated as effectively stable and **never decays → an empty histogram / zero deposits**.
  This bites *any* long-lived ion source, so it is easy to forget and hard to diagnose (the run
  "succeeds" with no counts).
- A **metastable isomer with a minutes-scale half-life is still not fully decayed** under these
  settings — e.g. Cs-137 → **Ba-137m (2.55 min)** → 662 γ emits in only ~half the primaries, so the
  662 peak is ~2× low. Use a parent that decays directly (Ba-133) as the reliable probe, or handle the
  isomer explicitly.
- Match the **world material** to what the MC assumes: for the electron-bearing channels (IC/β
  electrons crossing the source→detector gap) use an **air world** (`write_gdml(..., vacuum_world=false)`)
  so G4's electrons air-degrade like the engine's; a vacuum world lets them through undegraded. Photon
  air attenuation over cm-scale gaps is negligible either way.

**Harness flags** (`ceelo_g4val <gdml> <mac> <out.csv> [flags]`):
- `--cone-bias` — cone importance sampling (point sources only; matches the MC cone estimator).
- `--histogram` — write `<out>_histogram.csv` pulse-height spectrum (`# Events:` header; lower-channel-energy bins).
- `--correlated-gamma` — enable `G4PolarizationTransition` γ-γ angular correlation W(θ) (cascade summing).
- `--diagnostics` — 11 per-interaction crystal-side histograms (debug scatter angles / interaction depths / escape energies).
- `--entry-diag` — external-electron-entry counters at the crystal (the source-electron escape channel; u/s source-class decomposition printed when s-events exist).
- `--lowcut` — production cut 0.001 mm (default 0.7 mm); needed for fluorescence X-rays from heavy attenuators (e.g. Pb Kα range < 0.7 mm).
- `--em-standard` — use `G4EmStandardPhysics` instead of option4.
- `--fep-precision <p>` — run until FEP relative precision ≤ p (cap 100M events).

Physics: **FTFP_BERT + `G4EmStandardPhysics_option4`** (full Doppler broadening). Deexcitation
`SetFluo(true)`, `SetAuger(true)`, `SetDeexcitationIgnoreCut(true)` — the last is critical so
fluorescence X-rays from heavy-element attenuators are produced (else killed by the default cut).

**Inspecting G4 internals (how G4 does things):**
- **Source tree** `<geant4-source>/source/processes/electromagnetic/` — the EM models CeeLo is
  validated against: `G4LowEPComptonModel` (our Compton reference), `G4eCoulombScatteringModel` /
  `G4WentzelVIModel` / `G4UrbanMscModel` (multiple scattering — relevant to the source-escape tails),
  `G4eBremsstrahlung` / `G4SeltzerBergerModel`, `G4eplusAnnihilation`, `G4PEEffectFluoModel`,
  `G4Livermore*` (low-energy). Read these to compare a sampling/cross-section choice against ours.
- **Aggregate reference export:** `tools/geant4_validation/` records configured
  spectra, rates, efficiencies, and histograms for comparison. Raw G4
  cross-section columns are not production inputs and must not be committed.
  Diagnostic binaries may still print transient spot checks during model
  investigations (for example, the cfg-7 Pb K-edge investigation).
- **In-macro introspection:** add `/run/verbose 1`, `/process/em/printParameters`, `/process/list`,
  `/process/em/verbose 1` to a macro to print the active processes, models, and EM parameters G4 uses.

### Validated Configurations

| Config | Detector | Attenuator | Source | FEP agreement | Total agreement |
|--------|----------|------------|--------|---------------|-----------------|
| 1 | 3"x3" NaI | bare | point, 10cm on-axis | ≤ 0.9% | ≤ 0.6% |
| 2 | 3"x3" NaI | 1mm Al | point, 10cm on-axis | ≤ 0.7% (≤1173); −2.1% @ 3000 | ≤ 0.6% |
| 3 | 2"x2" LaBr3 | 0.5mm Al | point, 5cm on-axis | ≤ 0.6% | ≤ 0.5% |
| 5 | 1x1x0.5cm CZT | bare | point, 5cm on-axis | ≤ 1.3% (30-662 keV); −5.3% @ 1500 | ≤ 0.5% |
| 6 | 3"x3" NaI | bare | point, 15cm 45deg off-axis | ≤ 0.6% | ≤ 0.6% |
| 7 | 3"x3" NaI | 1mm Al + 2mm Pb | point, 15cm on-axis | ≤ 0.8% (≥ 200 keV) | ≤ 0.9% (≥ 200 keV)*** |
| 8 | 3"x3" NaI | 0.5mm Al | Marinelli beaker, water | ≤ 1.2% (≥ 100 keV); −2.3% @ 59* | ≤ 1.1%* |
| 11 | 3"x3" NaI | bare | point, 10cm + 0.5cm Fe shield | ≤ 0.6% (200-2000); −1.7% @ 3000** | ≤ 0.6%** |
| 12 | 3"x3" NaI | bare | 10x15x20cm SS304 box, cellulose | ≤ 0.2% (≥ 200 keV); −3.8% @ 59 | ≤ 1.4% (≥ 200 keV); −3.0% @ 59 |
| 25 | GEM35-70 HPGe coax, sharp edge | bare | point, 5cm on-axis | ≤ 0.34% | ≤ 0.10% |
| 26 | GEM35-70 HPGe coax, **bulletized** + round-tipped bore | bare | point, 5cm on-axis | ≤ 0.74% | ≤ 0.14%† |

**Measured Aug 2026** against the committed GEANT4 references, from
`tests/data/ceelo_reference/` regenerated at ~0.3% precision on the EPICS2023
photon data. A few tenths of a percent of every entry is Monte Carlo statistics.

**Configs 25/26 are a matched pair** (added Aug 2026 with bulletization support):
the same GEM35-70 HPGe coax with a sharp and a bulletized front edge, 4M
cone-biased G4 events per energy (12M for config 26 at 45 and 88 keV) and
CeeLo at ~0.03% precision below 122 keV, ~0.1–0.35% above, over 45–1332 keV. Each config agrees with GEANT4 on its own (above), but the
quantity the pair exists to measure is the **bulletization effect itself**,
Δ = 1 − ε(26)/ε(25), computed independently in each code:

| E (keV) | Δ FEP CeeLo | Δ FEP G4 | difference | z | Δ total CeeLo | Δ total G4 | difference | z |
|---------|-------------|----------|------------|---|---------------|------------|------------|---|
| 45 | 16.43% | 16.39% | +0.04 pp | 1.09 | 16.16% | 16.17% | −0.01 pp | −0.19 |
| 59.5 | 15.27% | 15.19% | +0.08 pp | 1.84 | 15.11% | 15.06% | +0.04 pp | 1.03 |
| 88 | 12.25% | 12.25% | −0.00 pp | −0.01 | 12.26% | 12.28% | −0.02 pp | −0.41 |
| 122 | 9.03% | 9.09% | −0.05 pp | −0.90 | 9.39% | 9.42% | −0.04 pp | −0.70 |
| 344 | 4.16% | 4.09% | +0.07 pp | 0.26 | 5.26% | 5.25% | +0.01 pp | 0.07 |
| 662 | 2.94% | 3.15% | −0.21 pp | −0.50 | 4.24% | 4.42% | −0.18 pp | −1.19 |
| 1332 | 2.03% | 2.47% | −0.44 pp | −0.81 | 3.86% | 3.84% | +0.02 pp | 0.14 |

† **The apparent config-26 low-energy residual was statistics plus a shared
offset — not bulletization.** The first pass showed config 26 running +0.16 to
+0.28% above GEANT4 at 45 and 88 keV where config 25 did not, which looked
geometry-specific. Re-running CeeLo at ~0.03% precision (two independent
seeds, 45–122 keV; the committed rows are those runs combined) resolved it:

- At **45 keV** the excess vanished — it was CeeLo statistics. The original
  run carried ~0.10% uncertainty and the "excess" was ~2.4σ of it. Config 26
  now sits −0.08% FEP / −0.01% total.
- At **88 keV** a real +0.135% FEP offset survives and sharpens (z = 3.5–3.9),
  but it is **identical in both configs** — cfg 25 +0.135%, cfg 26 +0.136%, a
  difference of +0.000 pp. It has nothing to do with the fillet. 122 keV
  carries a similar shared offset of about −0.1%.

So what remains is a ~0.1% CeeLo-vs-GEANT4 offset on germanium with energy
structure across 60–122 keV, common to sharp and bulletized crystals alike and
therefore outside the scope of this feature. The leading candidate is the
photon-evaluation difference already documented above: CeeLo is on EPICS2023
while GEANT4 11.4.0 ships EPICS2017, which produced exactly this kind of
energy-structured sub-percent difference for iron in config 12. Confirming it
for germanium means comparing µ/ρ(Ge) against NIST XCOM the way commit 1ad196c
did for Fe/Cr/Ni. Tracked in TODO.md; it is well inside the ~1% goal.

**What the pair does and does not establish.** GEANT4's transport is
independent, but the *solid* it navigates comes from `write_gdml()` in the same
change, built from the same arc parameterisation the tracer uses. So the pair
validates "CeeLo's tracer agrees with CeeLo's exported solid, under GEANT4
physics" — it does not, by itself, establish that either matches ANGLE's
`bulletizingRadius`. That last link rests on the Pappus closed form and the unit
tests in `tests/test_geometry.cpp` (chords against an independent
point-in-solid oracle) and `tests/test_geant4_export.cpp` (exported profile
against the traced solid).

Every point agrees within |z| ≤ 1.9, on an effect that runs from 16% down to
2%. The effect is far larger than the 1.28% of crystal volume the fillet
removes, and strongly energy dependent, for a geometric reason worth recording:
from this source position bulletization removes **17% of the rays that hit the
crystal at all** but only **2.4% of the total path length**, because the rays it
eliminates are short grazing corner clips. At 45 keV, where photons stop within
the first millimetres, losing those rays (and pushing the entrance surface
deeper for the rest) costs ~16%. By 1332 keV a 0.5 cm chord interacts only
~14% of the time against ~64% for an average chord, so dropping them costs only
a few percent. That, not the volume fraction, sets the high-energy limit.

**EPICS2023 vs GEANT4's EPICS2017.** CeeLo's photon magnitudes now come from a
*newer* evaluation than the one GEANT4 11.4.0 ships, so exact agreement is no
longer the expected outcome. Migrating from the previous Geant4-fitted
representation moved the mean |deviation from GEANT4| over the 174 gated
(config, energy, metric) points from 0.366% to 0.432%, with 70/174 points moving
closer and the rest further. The two largest movements are **config 12 @ 59 keV**
(FEP −0.05% → −3.76%, total +0.92% → −3.02%) and **config 11 @ 3000 keV**
(FEP −0.13% → −1.65%). The accepted migration gate was a 0.7173% p95 /
1.6447% worst contribution-weighted attenuation change, with independent NIST
XCOM spot checks favouring the direct EPICS values over the historical fit; see
`tools/prepare_cross_sections/reports/migration_acceptance.json`.

\* **Config 8 reference was regenerated at 32M on the CURRENT geometry (June 26 2026).** The previous
March `g4_cfg8_*_v2` reference was **stale**: it predated the June source-effect/biasing/cascade
toolchain (GDML export rebuilt Jun 11, harness Jun 25) and sat **+0.9% (662 keV total) to +6.1%
(59 keV total) ABOVE both current MC and a fresh G4 run on the current geometry**. Against the fresh
reference, MC agreed to ≤0.7% total and ≤1.1% FEP across 100–2614 keV (z ≤ 2.4) on the
then-current photon data; on EPICS2023 (Aug 2026) that is **≤1.2% FEP and ≤1.1% total**
over the same range, with 59 keV FEP at −2.30%;
`compare_validation.py` gates **all six energies** (59/100/200/662/1173/2614; FEP tol 2.5%, total 1.0%;
**no SKIPs**). The earlier "low-energy Marinelli transport deficit" (FEP −4 to −7%, total −2 to −5.6%)
and the "−2.5% albedo-recovery FEP deficit @662" were **largely stale-reference artifacts** — fresh
current-MC-vs-current-G4 is −0.17% FEP @662, −0.25% total @59. The one **genuine residual is FEP −2.1%
at 59 keV** (z=5): the G4-harness u/s decomposition shows the unscattered FEP `eps_u` matches G4 to
0.01–0.3% at every energy (geometry + self-attenuation are right; water μ = XCOM to 0.6%), but the
small *scattered* FEP stream `eps_s` is ~15% low — MC under-places forward Rayleigh / small-angle-Compton
photons in the ±1.5 keV window (shrinks to ≈0 by 200 keV; same family as the MC-delta-vs-G4-0–1-keV-smear).
A small real albedo-recovery effect (G4 recovers crystal-escaped e⁻/γ backscattering off the surrounding
water; MC's ≤3-bounce photon-only re-entry + one-way electron escape miss it) may persist but is ~10×
smaller than the previously documented −2.5%. See the reference CSV header for the regeneration record. Config 8's Marinelli electron-escape path
uses the regime-aware `walk_in_source_geometry` (water is light-Z → fully analog).

\*\* Config 11 GDML uses `rmin=1e-4 cm` for the Fe source shield sphere to avoid G4 navigation issues when a GPS source is at the exact center of a solid sphere (`rmin=0`). Configs 11/12 cover 200–3000 / 59–2614 keV; total agreement above the 1022 keV PP threshold requires the source-shield secondary channels (annihilation gammas, shield bremsstrahlung, electron escape — see Source Geometry section), which together contribute up to ~6% of total efficiency at 3 MeV. References upgraded June 2026: cfg 11 {2000, 3000} and cfg 12 {662, 1173, 2614} rows are 64M-event G4 runs (old-vs-new G4 z ≤ 1.3); the previously documented "−1.3% FEP at 3000" was an old-reference fluctuation (vs 64M: −0.4 to −1.1% across repeat runs, part of the high-energy SB/CSDA family). **Config 12 totals agreed to −0.5 to +0.8% across 59–2614 keV (z ≤ 1.5)** on the then-current photon data; on EPICS2023 (Aug 2026) 200–2614 keV is +0.01 to +1.32% while 59 keV moved to −3.02% (see the Validated Configurations note) after the June 2026 skin-escape-transmission fix (`ElectronCsda::walk_in_source_geometry`). The pre-fix excess was +1.8–2.5% ABOVE G4 at 662–1173 keV (z = 4–7) — root-caused (June 2026 residual triage) to the source-electron *escape* channel delivering ~5× G4's electron-arrival rate (MC 8.7e-5/event vs the G4 electron-entry-with-deposit rate 1.96e-5 ± 0.11e-5) for sub-MeV skin escape from dense shields (the Gaussian-core Molière walk lacks large-angle tails/backscatter). The fix is an analog surface-emergence (1 − albedo) survival gate, A(E, Z_eff), applied to the escape outcome only (bremsstrahlung is kept): it suppresses the near-surface τ≈0 sub-MeV escapers G4 backscatters, while saturating to ≈1 above ~1.5 MeV so the MeV channel is preserved. Post-fix electron channel (`TOT` − `--no-source-electrons`): cfg 12 @662 = 2.11e-5 ± 0.35e-5/event vs the G4 electron-entry target 2.0e-5 ± 0.1e-5 (z = 0.3, was 8.7e-5 pre-fix); cfg 11 @3000 = 3.03% of total vs the G4-implied +3.05% (×1.21 → ×1.0 vs G4). The channel must stay enabled (without it, cfg 11 @3000 total is −3.05% ± 0.18%, z = 17). See Known Limitations → "Source-electron skin escape".

\*\*\* Config 7 at 100 keV: ~−25–31% total discrepancy, excluded from the automated gate (`SKIP` in `profiling/compare_validation.py`). **Localized June 2026** to a G4-side effect: G4's *in-tracking* μ_total(Pb) is smeared across the K-edge (88 keV) relative to its own G4EmCalculator export, while the exported table matches NIST XCOM to 0.23% (5.536 vs 5.549 cm²/g at 100 keV). Measured via cone-biased uncollided-FEP ratios through the 2 mm Pb stack (Δμ/μ = ln(FEP_G4/FEP_MC)/τ): **+0.35% ± 0.15% at 80 keV (below edge, μ pulled up), −18% ± 1.5% at 90 keV (just above, pulled down), −4.4% ± 0.5% at 100 keV (recovering)** — the classic λ-table edge-interpolation signature; 73.5/75 keV agree to ±0.2–1%. Consequences: G4 transmits/escapes far more just-above-edge flux; the K-escape channel gains ~+10–15% from the deeper PE-production profile (forward K escape ∝ μ_PE/(μ₁−μ₂), front-face dominated — NOT e^(−μ₁t)). The Kα-specific ratio (G4/MC ×1.50, Kβ ×1.22) is additionally driven by the escaped Kα/Kβ ratio: G4's 0.82 matches a simple emission×attenuation estimate (BR_Kβ/Kα ≈ 0.30, Δμ₂t ≈ 1.4), MC's 0.67 is low — follow-up: verify MC's Pb K-line branching ratios/energies from EADL. G4's run-to-run scatter resolved: three runs {1.366 ± 0.006, 1.476 ± 0.025, 1.377 ± 0.013}e-5; the high one was the outlier; G4 ≈ 1.371e-5 ± 0.6%. NOT an MC transport bug (MC tables = XCOM; below-edge transmission matches G4). Exclusion retained.

G4 reference data: 64M events/energy (cone-biased for point sources), 16M events/energy (isotropic for extended/shielded sources, except cfg 11 {2000, 3000} and cfg 12 {662, 1173, 2614} which are 64M — June 2026). MC comparison: 4M events/energy (or precision-targeted ~0.3% with auto biasing).

### Spherical / annular source geometries (June 2026, unnumbered)

Added with the cascade geometry-parity work (`examples/sphere_g4_validation.cpp`,
`build/examples/sphere_g4val/`). Three 3"×3" NaI geometries — **G-A** bare self-attenuating Thorium
sphere (R=2 cm), **G-B** trace source in a soil sphere (R=3 cm) + Al(0.2)/Fe(0.5) shells, **G-C**
void-center soil shell [2,3] cm + Fe(0.5) — at near (~2 cm gap) and far (~50 cm); Th-232 chain energies.
Vacuum world. G4 = 4M-event isotropic GPS volume source.

- **FEP agrees across all geometries/energies**: |z| ≤ 1.4 at 238.6/583/911 keV (e.g. G-A FEP z = −0.5/0.4/0.7; G-B/G-C within ±1.5). FEP is the priority metric and validates the geometry, sampling, and self-attenuation path. (G-A FEP at 2614 keV is −4.2%, z≈−3.3 — the high-energy SB/pair/CSDA family, amplified by the thick high-Z source, not geometry.)
- **Low-Z trace-source totals agree**: G-B and G-C total efficiency match G4 to ≲1.5% (|z| mostly ≤2) across 238.6–2614 keV — validates the trace + void-center-shell geometry and shielding.
- **High-Z self-attenuating total is low**: G-A (Thorium) total is **−22% / −4.1% / −2.6%** at 238.6/583/911 keV — the source-material **K-fluorescence** that is not emitted (Th Kα ≈90/93 keV escaping the source), confirmed by the G4 histogram and FEP-clean. See Known Limitations. (G-A 2614 keV total −4.0% is the coherent FEP+total high-energy family, separate from fluorescence.)

The far (~50 cm) MC efficiencies converge well via the auto two-stream/cone direct stream; a matching G4 cross-check at 50 cm needs much higher isotropic stats (cone bias is invalid with source scatter) and is left for a dedicated high-stats run. The cascade path on a sphere is exercised by `tests/test_spherical_source.cpp` (`co60_summing_out_on_sphere`); `compute_cascade` is geometry-agnostic (routes through `sample_source_position` + `transport_source_photon`).

### Generating Reference Data

```bash
cd build/examples
./generate_all_spectra [n_mc_histo] [n_g4_macro] [n_mc_multi]
# Defaults: 4M histogram, 100M G4 macro, 100k efficiency

# Then run G4 (generated script handles cone bias per config):
bash run_g4_all.sh
```

### Diagnostic Tools

- `diag_interactions [energy_keV] [num_events]` -- per-interaction MC histograms (11 CSVs)
- `ceelo_g4val ... --diagnostics` -- per-interaction G4 histograms
- Useful for debugging scatter angle distributions, interaction depths, escape energies

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

- Use lower channel energy (not bin center)
- `RealTime` line has no `#` prefix (compatible with spectral analysis software)
- `RealTime = N_events / (1E6 * omega_frac)` for cone-biased; `N_events / 1E6` for isotropic
- This represents equivalent measurement time for a 1 MBq source with 100% branching ratio

## Known Limitations

Open limitations, shortcomings, and deferred features are tracked in **`TODO.md`** (repo root) — the
single source of truth. Add new shortcomings there, not here. Detailed root-cause and fixed-issue history
is in the git history.

Standing notes that are easy to trip over:
- **High-energy FEP family.** Above ~2 MeV the FEP residual is negative and *coherent in sign across
  NaI/LaBr3 configurations*, growing with energy. Measured on the EPICS2023 data (Aug 2026): −2.06%
  (cfg 2 @3000), −1.65% (cfg 11 @3000), −0.62% (cfg 1 @2000), −0.21% (cfg 3 @3000), −0.06%
  (cfg 12 @2614). **Below ~1.3 MeV the residual is small** (|Δ| ≤ 0.6% for cfgs 1/2/3/6/11/12 at
  0.3–0.5% combined σ). The coherence across configurations points at the Seltzer-Berger
  bremsstrahlung + CSDA electron-escape approximations rather than anything configuration-specific.
  This supersedes the older, vaguer "−1 to −2% FEP residual, root cause unknown". CZT (cfg 5) is
  separate: −2.28/−4.12/−5.34% at 800/1000/1500 keV from thin-crystal electron escape.
- **Config 12 at 59 keV is a data difference, not a transport defect** (attributed Aug 2026;
  probe in the dev-only `studies/xs_probe/`). cfg 12 @59 keV moved from −0.05% to −3.76% FEP
  (−3.72% ± 0.42%, z = 8.9) across the EPICS2023 migration. Measuring μ at 59 keV on both data
  sets: **SS304 +1.730%**, cellulose +0.249%, NaI −0.014%. The 0.2 cm steel wall has
  τ = μ·t = 1.99, so first-order attenuation predicts
  ΔFEP = −(τ·sec θ·Δμ_SS + τ_cell·Δμ_cell) = **−3.79%**, against **−3.72% ± 0.42%** observed —
  agreement to 0.2σ. The same model reproduces the energy dependence (predicted/observed:
  −0.24%/−0.74% at 100 keV, −0.12%/+0.00% at 200 keV, ≈0 above), i.e. the effect is confined to
  the steel's photoelectric regime.

  The underlying change is a **smooth ~30 keV-wide bump in the iron photoelectric cross section**,
  peaking at **+2.23% near 58 keV** and falling to +0.34% at 40 keV and +0.81% at 90 keV. It is
  smooth rather than oscillatory, so it is an evaluation-level difference between EPICS2023 and the
  previous Geant4-fitted representation, not a grid/interpolation artifact. 59 keV sits almost
  exactly on the peak — the configuration is unluckily placed.

  **Arbitrated against NIST (Aug 2026): EPICS2023 is the correct one.** Against the NIST *X-Ray
  Mass Attenuation Coefficients* tables (Hubbell & Seltzer, SRD 126 — derived from neither
  EPICS2017 nor EPICS2023), mu/rho deviations are:

  | element | 50 keV | 60 keV | 80 keV | 100 keV |
  |---|---|---|---|---|
  | Fe pre-migration | −1.30% | −1.26% | −0.36% | +0.08% |
  | Fe EPICS2023 | **+0.31%** | **+0.36%** | **+0.39%** | **+0.32%** |
  | Cr pre-migration | −1.67% | −1.32% | −0.34% | +0.11% |
  | Cr EPICS2023 | **+0.31%** | **+0.35%** | **+0.33%** | **+0.27%** |
  | Ni pre-migration | −0.97% | −1.07% | −0.42% | +0.06% |
  | Ni EPICS2023 | **+0.31%** | **+0.38%** | **+0.35%** | **+0.29%** |

  At 50–60 keV EPICS2023 is ~4x closer to NIST. Its offset is also *flat* (+0.27 to +0.39% over
  3 elements x 4 energies), the signature of a coherent evaluation, whereas the old data carries
  strong energy structure (−1.3% at 50–60 keV recovering to ~0 by 100 keV) — the same ~30 keV bump
  seen in the direct comparison. The ~1.6 percentage-point mu difference times tau ~ 2 predicts
  ~3.2% in transmitted flux, consistent with the −3.76% observed.

  **So config 12 @59 keV is CeeLo being more accurate than GEANT4, not less.** Full working in the
  dev-only `studies/xs_probe/nist_arbitration.md`.

- **Config 8 at 59 keV is unchanged and separate.** Its shift across the migration is −1.14% ±
  0.59% (z = 1.9) — statistically marginal, and consistent with the ≈−2% residual documented for
  that point since June 2026. Water's μ moved only +0.264% at 59 keV and the beaker's water path
  gives τ ≲ 1, so attenuation accounts for ≲0.3% of it. This remains the scattered-FEP stream
  `eps_s` (small-angle Compton and forward Rayleigh that should land in the ±1.5 keV window)
  under-populated relative to GEANT4.
- **Doppler broadening** (impulse approximation, default ON) broadens the Compton *edge* to match
  GEANT4; the subshell binding energy is deposited locally, so FEP/total efficiency are statistically
  unchanged. **Never cite Doppler broadening as a cause of FEP/total discrepancies** — its effect on
  those is well below the FEP window.
- **Source-electron skin escape** uses an empirical surface-emergence albedo gate
  (`ElectronCsda::source_escape_survival_exit`), now **regime-aware**: evaluated at the EXIT state
  (exit KE + exit-material Z), with a light-Z floor (η₀ ≤ 0.16 → fully analog, so Al positrons / water
  match G4 with no per-caller flag) and an exit-energy window. This one material-general gate replaced
  the old birth-energy gate + the `apply_skin_escape_gate=false` positron flag + the Fe-only `kAlbScale`.
  Validated across Z = 5→82 (PE/water/Al/Fe/SS304/Sn/W/Pb) to ~1–2% vs GEANT4. The MeV electron-entry
  channel it gates is REQUIRED for cfg 11/12 totals above the 1022 keV PP threshold — disabling it gives
  −3% at 3 MeV.
  - The escape model is **COMPILE-TIME selectable** (`cmake -DCEELO_SOURCE_ESCAPE_MODEL=gate|tails|gs`,
    default `gate`): `gate` = this empirical gate; `tails` / `gs` are first-principles (straggling +
    screened-Rutherford hard tail, no gate) differing only in the soft Gaussian-core width
    (Highland-minus-hard-variance vs the screened-Rutherford soft transport moment). They differ only
    inside `walk_in_source_geometry` (per-step MSC + the escape decision):

    | Model | Gaussian core θ₀ | hard tail | straggling | escape decision |
    |-------|------------------|-----------|------------|-----------------|
    | **gate** (0) | Highland (per-step) | — | — | empirical exit-state (1−albedo) gate |
    | **tails** (1) | Highland-B2 − hard variance | screened-Rutherford μ<μ_c | Bohr/clamped | analog (walk decides) |
    | **gs** (2) | √(G₁,soft·Δx), soft transport moment of the screened cross-section | screened-Rutherford μ<μ_c | Bohr/clamped | analog (walk decides) |

    Total-efficiency error vs GEANT4 at 662 keV (%, MC ~0.2%), across the Z = 5→82 box/shell sweep:

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

    **Verdict (mean |err| @662): gate 0.56% < tails 0.62% < gs 0.92%** — `gate` is the validated default.
    `tails` is competitive and better at extreme Z (the gate's 0.92 albedo cap over-suppresses Pb), but its
    per-step Gaussian+tail over-diffuses sub-MeV mid-Z **boxes** (+0.6–1.2% vs gate; the double-counting
    variance subtraction halved that) — a **box-specific multi-wall artifact**, since a large-angle scatter
    is as likely to redirect toward a wall as away. On sphere/shell geometries (G-B/G-C at 238.6–2614 keV,
    and the cfg-11 Fe shell) all three models agree to **~±1%**, so the first-principles models are as good
    as the gate for radial escape; the gate stays default because it is the only one that also handles the
    multi-wall box case well. `gs` under-scatters ~30% (its screened-Rutherford G₁ omits the Mott spin and
    nuclear-form-factor corrections GEANT4 uses), so its soft core is too narrow and boxes over-escape more;
    a proper Mott / tabulated Goudsmit-Saunderson MSC is the path to retiring the gate (`TODO.md` B2/B3).
    **Timing**: `tails` ≈1.03–1.10× `gate`, `gs` ≈1.1–1.2× — immaterial, since the source-escape walk is
    <10% of runtime (photon transport dominates and is identical across models).
    High-Z is **not** a single gate knob: Sn/W over-predict while Pb under-predicts at the *same* capped
    albedo, and lowering the cap 0.92→0.85 worsened Sn (+0.32→+0.84) and W without helping Pb (whose
    electron channel is tiny — Pb stops electrons — so its −1.46% is mostly photon transport).
- **Source-material fluorescence is not emitted**: photoelectric absorption in the source material/shields
  drops the photon without a characteristic K X-ray, so **total** efficiency for **high-Z
  self-attenuating** sources is low where PE dominates (validated June 2026: a Thorium sphere is −22% /
  −4% / −2.6% in total at 238.6 / 583 / 911 keV vs GEANT4, tracking the missing Th Kα ≈90/93 keV escape).
  **FEP is unaffected** (the X-ray cannot fall in the primary peak) and **low-Z sources are unaffected**.
  Tracked in `TODO.md` (Photon transport physics → "Source-material fluorescence not emitted").

## Common Pitfalls

- **Physics tables are auto-generated**: Do not edit `element_data.cpp`,
  `photon_epics_data.cpp`, `relaxation_epics_data.cpp`, or
  `estar_stopping_data.cpp` by hand. Use the locked generators and update the
  generated manifest last.
- **C++ `const` linkage**: Namespace-scope `const` arrays have internal linkage by default. If a header declares `extern const T arr[]`, the .cpp definition must also use `extern`.
- **Eigen `.cross()`**: Requires `#include <Eigen/Geometry>`, not just `<Eigen/Core>`.
- **Batch size in `compute()`**: Must be much smaller than `num_events` so threads check termination frequently. Currently `max(10000, num_events/10)`.
- **Importance sampling denominator**: The IS estimator is `(1/N) * sum(w_i * I_i)`. Using `sum(w_i)` in the denominator gives intrinsic efficiency (cancels the cone weight), not absolute efficiency.
- **Source-shield GDML coincident surfaces**: For an extended source (cylinder/box) wrapped by a shield, the shell's subtracted inner solid must be inflated by a tiny epsilon (`Geant4Export.cpp`, `kZeroDimEps = 1e-4 cm`) so its inner face does not coincide with the source-material outer face (or the previous shell). Coincident G4 boolean surfaces stall navigation ("track stuck", ~20k warnings/run) and — critically — **spuriously inflate the crystal electron-entry `--entry-diag` rate ~1.75×** while leaving totals (geometry-robust) unchanged, so they silently corrupt channel-level comparisons. Fixed June 2026; `benchmark_mc_configs` also now exports source-effect box configs (12/20/21/22) with `vacuum_world=true` to match `generate_all_spectra`. A residual ~1.5k warnings/run at `active_crystal_PV` is a separate benign detector artifact present in the validated point-source configs too.
