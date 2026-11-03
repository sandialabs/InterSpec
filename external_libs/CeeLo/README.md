# CeeLo

**The everyday spectroscopist's Monte Carlo efficiency estimator.**

CeeLo is a small, self-contained C++17 Monte Carlo photon-transport library that
answers one question well: *given this detector, this shielding, and this source,
what fraction of emitted photons will I actually count?*

Give it a detector geometry (NaI, LaBr₃, HPGe, CZT, …), any attenuators or
shielding around it, and a source configuration, and it computes:

- **Full-energy peak (FEP) efficiency** — probability a photon deposits 100 % of
  its energy in the crystal
- **Total efficiency** — probability it deposits any energy at all
- **Pulse-height spectrum** — energy-deposit histogram, in a format spectral
  analysis software can read. `enable_spectrum_fidelity(true)` adds the detector physics
  that shapes the spectrum without moving FEP or total: x-rays from Compton K-shell
  vacancies, which fill out the K x-ray escape peaks above ~300 keV
- **True-coincidence (cascade) summing corrections** — summing-in and
  summing-out factors for multi-gamma decays

Sources can be points, cylinders (solid or annular), boxes (solid or hollow
shells), spheres (solid or void-centre shells), or Marinelli beakers — each with
its own material, self-attenuation, concentric shielding, and optional
attenuating cores inside hollow shapes.

---

## ⚠️ Status: under active development

**This library is not finished, and you should not rely on it.**

- Interfaces and numerical results may change without notice.
- It has been validated only over the configurations and energy range in the
  table below. Outside those, its accuracy is unknown.
- Known physics gaps are real and documented — see
  [Shortcomings](#shortcomings) and [TODO.md](TODO.md).
- **Do not use it for regulatory, safety, dosimetric, or any other
  decision-critical work.**

---

## Why it exists

InterSpec needed a Monte Carlo based method to help calculate total detector
efficiency for cascade-decay corrections. Bringing in GEANT4 or MCNP was overkill
and would introduce distribution difficulties, and requiring users to run their
own simulations is unreasonable. GEANT4 and MCNP are powerful general purpose
simulation toolkits, but for gamma spectroscopy only a very small subset of their
functionality is needed.

- A GEANT4 installation is **48 MB of libraries plus 2.1 GB of data**
  (726 MB of it just the low-energy EM data set).
- General-purpose transport cannot assume anything about your problem, so it
  cannot apply the variance reduction that a detector-efficiency calculation
  invites — the source is small, the detector subtends a known solid angle, and
  the quantity of interest is a peak area.

CeeLo is the narrow-purpose alternative: **a ~1.9 MB static library (Release)
with every cross-section table compiled in and no runtime data files at all**,
built around a biasing policy written specifically for this problem. It is fast
enough to run interactively inside a spectroscopy tool, which is where it is used.

The trade is deliberate: CeeLo does photon transport (with the electron physics
needed to get photon deposits right) and nothing else — no neutrons, no hadronic
physics, no optical photons, no beamline. It is validated from 30 keV to 3 MeV;
its photon data span 10 keV – 20 MeV.

---

## Relationship to InterSpec

CeeLo is developed alongside [InterSpec](https://github.com/sandialabs/InterSpec),
Sandia's gamma-spectroscopy analysis application, and is vendored into it.

**It has no dependency on InterSpec** and is usable on its own. The library
needs only a C++17 compiler, CMake, and Eigen3; Boost is needed only to build
the tests. Everything InterSpec-specific lives on the other side of the boundary.

---

## Build

```bash
cmake -B build -S .
cmake --build build -j8
cd build && ctest --output-on-failure
```

The default build type is `RelWithDebInfo`. The default suite registers ~32 ctest
targets (up to ~43 with the optional features below) and takes ~15–30 minutes:
these are Monte Carlo simulations with statistical convergence checks, not unit
tests of arithmetic.

**Dependencies:** C++17 compiler, CMake ≥ 3.15, Eigen3, and Boost (Boost.Test,
for the tests only).

**CMake options:**

| Option | Default | What it does |
|---|---|---|
| `CEELO_BUILD_TESTS` / `CEELO_BUILD_EXAMPLES` | ON | Boost.Test suites / runnable examples and benchmark drivers |
| `CEELO_NATIVE_ARCH` | ON standalone, OFF as a subproject | Adds `-march=native`; turn it off for binaries that must run on other CPUs |
| `CEELO_WITH_SANDIADECAY` | OFF | Cascade-summing adapter driven by real decay data ([SandiaDecay](https://github.com/sandialabs/SandiaDecay)); point `CEELO_SANDIADECAY_DIR` at its source tree |
| `CEELO_BUILD_RESPONSE` | OFF | Storable, MC-parameterized detector-response payload + EFFTRAN-style efficiency transfer; needs header-only RapidXml (deliberately not vendored; `CEELO_RAPIDXML_DIR`) |
| `CEELO_RUN_MC_TESTS` | OFF | Registers the slow MC-driven response tests (with `CEELO_BUILD_RESPONSE`) |
| `CEELO_BUILD_GEANT4_VALIDATION` | OFF | The GEANT4 comparison harness (needs GEANT4) |
| `CEELO_SOURCE_ESCAPE_MODEL` | `gate` | Source-electron escape model: `gate` \| `tails` \| `gs` |

## Quick start

```cpp
#include "efficiency/EfficiencyCalculator.h"
#include "materials/Material.h"
#include <iostream>

int main() {
    using namespace ceelo;

    Material nai = make_NaI();

    EfficiencyCalculator calc;
    // 3"x3" NaI(Tl): a radius (HALF the diameter) and the FULL crystal length.
    // See CRYSTAL DIMENSION CONVENTION in src/geometry/Geometry.h.
    calc.set_detector(&nai, CylinderDims{3.81, 7.62});
    calc.set_point_source(Eigen::Vector3d(0.0, 0.0, -10.0));  // 10 cm from the face

    auto r = calc.compute(661.7 /*keV*/, 100'000 /*events*/, 0 /*auto threads*/);

    std::cout << "FEP   = " << r.full_energy_peak_efficiency
              << " +/- " << r.fep_uncertainty << "\n"
              << "total = " << r.total_efficiency
              << " +/- " << r.total_uncertainty << "\n";
}
```

Every result carries its statistical uncertainty. See [`examples/`](examples/)
for extended sources, shielding, spectra, and cascade summing.

---

## How it is validated

CeeLo is validated against **GEANT4 11.4.0** (`FTFP_BERT` +
`G4EmStandardPhysics_option4`), and the comparison is set up so the two codes
cannot silently disagree about the problem: **the GDML geometry and the GPS
source macro that GEANT4 runs are exported by CeeLo itself**, so geometry and
source are identical by construction rather than by careful transcription.

The GEANT4 results are committed in
[`tests/data/geant4_reference/`](tests/data/geant4_reference/), so anyone can
re-run the comparison without installing GEANT4: regenerate CeeLo's side as
described in [`tests/data/ceelo_reference/README.md`](tests/data/ceelo_reference/README.md),
then run

```bash
python3 profiling/compare_validation.py
```

### Validated configurations

CeeLo against 16M–64M-event GEANT4 references, both at ~0.3 % Monte Carlo
precision, so a few tenths of a percent of every entry is statistical. Bounds are
the largest |deviation| over the gated energies (30 keV – 3 MeV unless noted).

| # | Detector | Attenuator | Source | FEP | Total |
|---|---|---|---|---|---|
| 1 | 3″×3″ NaI | bare | point, 10 cm on-axis | ≤ 0.7 % | ≤ 0.5 % |
| 2 | 3″×3″ NaI | 1 mm Al | point, 10 cm on-axis | ≤ 0.9 % | ≤ 0.5 % |
| 3 | 2″×2″ LaBr₃ | 0.5 mm Al | point, 5 cm on-axis | ≤ 0.8 % (≤ 1332 keV); +0.6 to +0.8 % at 2–3 MeV | ≤ 0.4 % (≥ 50 keV) |
| 5 | 1×1×0.5 cm CZT | bare | point, 5 cm on-axis | ≤ 0.8 % (≤ 800 keV); −1.5 % @ 1000, −2.5 % @ 1500 (\|z\| ≤ 1.7) | ≤ 0.6 % |
| 6 | 3″×3″ NaI | bare | point, 15 cm, 45° off-axis | ≤ 0.5 % | ≤ 0.8 % |
| 7 | 3″×3″ NaI | 1 mm Al + 2 mm Pb | point, 15 cm on-axis | ≤ 0.8 % (≥ 200 keV) | ≤ 0.5 % (≥ 200 keV) |
| 8 | 3″×3″ NaI | 0.5 mm Al | Marinelli beaker, water | ≤ 0.9 % (≥ 100 keV); −1.6 to −1.9 % @ 59 | ≤ 0.6 % |
| 11 | 3″×3″ NaI | bare | point, 10 cm + 0.5 cm Fe shield | ≤ 0.8 % (200–3000 keV) | ≤ 0.8 % |
| 12 | 3″×3″ NaI | bare | 10×15×20 cm SS304 box, cellulose | ≤ 0.2 % (≥ 200 keV); −4 % @ 59 | ≤ 1.4 % (≥ 200 keV); −3 % @ 59 |
| 25 | GEM35-70 HPGe coax, sharp edge | bare | point, 5 cm on-axis (45–1332 keV) | ≤ 0.9 % | ≤ 0.4 % |
| 26 | GEM35-70 HPGe coax, bulletized | bare | point, 5 cm on-axis (45–1332 keV) | ≤ 1.1 % | ≤ 0.4 % |
| 27/28 | GEM35-70 HPGe coax, bulletized | bare | point, 2 / 10 cm, 0.5 cm Fe shell | +0.1 to +0.4 % @ 122; −8 to −9 % @ 60 | ≤ 0.2 % @ 122 |

What lies beyond ~1 %:

- **CZT above 1 MeV** — the GEANT4 reference itself carries 1–1.4 % there; the
  deviations are within 1.7σ.
- **59–60 keV behind or inside iron** (configs 12 and 27/28) — a photon-*data*
  difference, not transport. CeeLo uses **EPICS2023**, GEANT4 11.4.0 ships
  **EPICS2017**, and the two differ by ~2 % in iron's photoelectric cross section
  near 58 keV. Against the independent NIST *X-Ray Mass Attenuation
  Coefficients* tables, **CeeLo is the closer of the two** (+0.3 % for Fe, Cr, Ni
  at 50–60 keV, against −1.3 %). Through 2 mm of steel (τ ≈ 2) that becomes ~4 %
  in transmitted flux; through 0.5 cm of iron at 60 keV, ~8 %. It vanishes by
  200 keV. These 60/88 keV rows are skipped in the automated gate.
- **The water Marinelli beaker at 59 keV** — small-angle Compton and forward
  Rayleigh that should land inside the peak window are slightly
  under-recovered.
- **Config 12 total, +0.9 to +1.3 % at 200–1173 keV** — part of a small,
  coherent extended-source total offset tracked in TODO.md.

One documented exclusion: **config 7 at 100 keV**, where a ~−25 to −31 % total
discrepancy was traced to a GEANT4-side artifact (its in-tracking μ_total for Pb
is smeared across the 88 keV K-edge relative to its own exported tables, which
match NIST XCOM to 0.23 %).

Spherical and annular sources, sources with attenuating cores, and the
cascade-summing path are validated separately — see [DESIGN.md](DESIGN.md).

---

## Performance compared to GEANT4

Both codes were run on the **same exported geometry and source**, both
cone-biased, both stopping when the **FEP efficiency reached 0.3 % relative
statistical precision** — so this is time-to-equal-answer, not time-per-event.
Single-threaded, Apple M1 Pro, Apple clang 21, `-O2 -march=native`, measured
Aug 2026.

CPU seconds, CeeLo / GEANT4, and the speed-up:

| Config | 59 keV | 662 keV | 1332 keV | 2614 keV |
|---|---|---|---|---|
| 1 — NaI 3″×3″, bare | 0.05 / 2.11 (**43×**) | 0.92 / 8.16 (**8.9×**) | 2.14 / 15.30 (**7.2×**) | 5.77 / 33.58 (**5.8×**) |
| 3 — LaBr₃ 2″×2″ + Al | 0.05 / 2.16 (**45×**) | 1.26 / 10.41 (**8.2×**) | 2.90 / 19.55 (**6.7×**) | 7.90 / 43.17 (**5.5×**) |
| 5 — CZT 1×1×0.5 cm | 0.07 / 3.48 (**49×**) | 2.59 / 23.18 (**9.0×**) | — | — |
| 7 — NaI + Al + 2 mm Pb | 1.25 / 81.78 (**66×**, at 200 keV) | 2.46 / 23.03 (**9.4×**) | 4.76 / 37.61 (**7.9×**) | 12.61 / 81.48 (**6.5×**) |

**Median speed-up 8.6×; 5.5–9.4× across most of the range, 43–66× at low energy
and behind shielding.** Subtracting GEANT4's ~1.0 s fixed initialization cost the
median is still 7.7×.

Two caveats:

- **Both codes needed essentially the same number of events** (median ratio
  1.03), so this is per-event cost, not a cleverer estimator — with one
  exception: for the thin CZT crystal at 662 keV, CeeLo's forced-first-interaction
  biasing got there in 0.9M events against GEANT4's 5.4M.
- **CeeLo is solving a much narrower problem.** Its speed comes partly from
  geometry-specific assumptions a general-purpose code is not allowed to make.
  This is not a claim that CeeLo is a better transport code than GEANT4; it is a
  claim that a specialised tool is cheaper for the one job it does.

For shielded and extended sources GEANT4 cannot use cone biasing at all (source
scatter can redirect photons from any initial direction), so its reference runs
needed 16–64M isotropic events per energy — a further cost multiplier the
point-source table above understates.

### Distribution size

| | CeeLo | GEANT4 11.4.0 |
|---|---|---|
| Library | `libCeeLo.a` — **~1.9 MB** (`Release`) | `lib/` — **48 MB** |
| Runtime data files | **none** | `G4EMLOW` — **726 MB** (EM only; full data set 2.1 GB) |
| Physics data, compiled in | **~610 KiB** of constant data | n/a — loaded from disk |
| Self-contained example binary | **~1.0 MB** stripped | — |

Library size is for a `Release` build without the optional features. The default
`RelWithDebInfo` archive carries debug symbols and is ~27 MB; those symbols do
not reach a linked, stripped binary.

CeeLo has no runtime data directory, no environment variables to set, and
nothing to install alongside the binary.

---

## Where the data comes from

All physics data is compiled into the library. Full citations, licences, and
regeneration instructions are in
[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).

| Data | Source |
|---|---|
| Photon cross-sections (photoelectric, K-shell photoelectric, Compton, Rayleigh, pair), Z=1–98 (H–Cf), 10 keV–20 MeV | Direct **EPICS2023 EPDL** records, with process-specific grids retaining thresholds and absorption-edge sides |
| Incoherent scattering functions S(x,Z), coherent Rayleigh sampling CDFs | Direct **EPICS2023 EPDL**; form factors are used during offline CDF generation but are not embedded as runtime lookup tables; Compton-profile support and atomic weights via BSD-licensed [xraylib](https://github.com/tschoonj/xraylib) 4.2.1 |
| Bremsstrahlung spectra (Seltzer–Berger), Z=1–92 | Public-domain **NIST EPQ** at pinned commit `ce1e589`, with one fixed low-κ interpolation rule; Np–Cf reuse the uranium table |
| Atomic relaxation: K and L fluorescence yields, lines, Coster–Kronig (Z=1–99 for decay daughters) | Direct **EPICS2023 EADL** MF=28/MT=533 |
| Electron stopping powers | Direct **NIST ESTAR** collision and radiative values for every Z=1–92 at all 53 nodes from 10 keV–20 MeV; Np–Cf reuse the uranium values |
| Material compositions | **ICRU Report 37** (air), **IAEA-TECDOC-1011**, CRC Handbook |
| Nuclear decay data (cascade summing, optional) | **ENSDF** via [SandiaDecay](https://github.com/sandialabs/SandiaDecay) |

The bremsstrahlung generator builds the GEANT4-aligned 32-node κ grid from NIST
EPQ by general interpolation of public-domain data (κ=0.025 and 0.075 from
low-order polynomials through the neighbouring EPQ nodes; every other node copied
as is) — not extraction or fitted reproduction of G4EMLOW coefficients.

---

## Shortcomings

Most important first. TODO.md has the full, prioritised list.

**Known accuracy gaps**
- **Only K-shell x-rays are emitted.** L lines are deposited locally, so L-escape
  peaks are missing and a bare self-attenuating thorium sphere still reads ~1 % low
  in total efficiency against GEANT4 (it read −22 % at 238.6 keV before source K
  fluorescence was added).
- **Air between the source and the detector is not transported by default.**
  An opt-in analytic attenuation of the unscattered photons exists
  (`set_air_attenuation(AirAttenuation::AnalyticNoScatter)`), but it is not yet
  validated against GEANT4, and air scatter is not modelled at all.
- **Source-electron skin escape** uses an empirical, GEANT4-calibrated
  exit-state albedo gate rather than a full Mott / Goudsmit–Saunderson multiple
  scattering treatment. It is validated across Z = 5–82 to 1–2 %, but it is a
  calibration, not first principles. Two first-principles alternatives are
  selectable at compile time (`CEELO_SOURCE_ESCAPE_MODEL=tails|gs`) and are
  currently slightly less accurate on box geometries.
- **High-energy FEP** runs a few tenths of a percent high at 2–3 MeV (LaBr₃
  +0.6 to +0.8 %, NaI up to +0.7 %; individual rows within ~2.5σ).
- **Thin crystals above 1 MeV**: CZT FEP reads up to −1.5 % / −2.5 % at
  1000 / 1500 keV, within the GEANT4 reference's own statistics.
- **Low-energy extended sources**: a water Marinelli beaker runs −1.6 to −1.9 %
  FEP at 59 keV; ≤ 0.9 % from 100 keV.
- **Efficiency transfer** (`CEELO_BUILD_RESPONSE`) is ≤ 1 % on-axis in the far
  field but a few percent near-field and off-axis; queries carry inflated
  uncertainties and provenance flags in those regimes.

**Scope**
- Validated roughly 30 keV – 3 MeV. The photon data cover H–Cf (Z = 1–98) from
  10 keV to 20 MeV, but above uranium the electron tables are uranium's (a
  few-percent approximation for electron energy loss, up to 9 % for Cf
  bremsstrahlung), and nothing above 3 MeV is in the GEANT4 gate. No photonuclear
  reactions, no neutrons, no optical photon transport, no hadronic physics.
- Only the geometries in the validated-configuration table have been checked
  against GEANT4. Everything else is extrapolation.

**Statistical honesty**
- Reported uncertainties are statistical only. They do not include
  cross-section data uncertainty, and they do not include the modelling
  approximations above. When CeeLo says ±0.3 %, that is the Monte Carlo
  precision, not the accuracy.

---

## Repository layout

```
src/            The library
  cascade/        True-coincidence summing (analytic + correlated MC)
  cross_sections/ Compiled-in photon data and its interpolation
  efficiency/     EfficiencyCalculator — the main entry point
  export/         GDML + GEANT4 macro export
  geometry/       Ray tracing, attenuator shells, source geometry
  io/             Storable detector response, efficiency transfer (mostly optional)
  materials/      Material definitions and built-ins
  physics/        Electron CSDA, stopping powers, multiple scattering
  transport/      Photon transport: Compton, PE, Rayleigh, pair, fluorescence
examples/       Runnable examples and benchmark configurations
tests/          Boost.Test suites + reference data
tools/          Cross-section regeneration, GEANT4 validation harness
profiling/      The GEANT4 validation gate script
visualization/  Dependency-free WebGL viewer for exported GDML geometries
```

Further reading:

- **[DESIGN.md](DESIGN.md)** — physics choices, biasing policy, geometry model,
  the GEANT4 validation procedure, and the standing explanations of known
  differences.
- **[TODO.md](TODO.md)** — the open-items register, in rough order of priority.
- **[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md)** — dependencies, data
  provenance, and the published methods implemented from the literature.
- **[visualization/README.md](visualization/README.md)** — how to eyeball an
  exported geometry before running a long simulation.

---

## Licence

CeeLo is distributed as part of InterSpec under the **GNU Lesser General Public License, version
2.1** — see [LICENSE.txt](LICENSE.txt).

Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

Questions, bug reports, and attribution concerns: wcjohns@sandia.gov or
interspec@sandia.gov.
