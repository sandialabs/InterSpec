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
  analysis software can read
- **True-coincidence (cascade) summing corrections** — summing-in and
  summing-out factors for multi-gamma decays

Sources can be points, cylinders (solid or annular), boxes (solid or hollow
shells), spheres (solid or void-centre shells), or Marinelli beakers — each with
its own material, self-attenuation, and concentric shielding.

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
efficiency for cascade-decay corrections.  Bringing in GEANT4 or MCNP was overkill
and would introduce distribution difficulties, and requiring users to run their
own simulations is unreasonable.  GEANT4 and MCNP are powerful general purpose
simulation toolkits, but for gamma spectroscopy, we only need a very small subset
of functionality.  This is a simple solution that can be integrated into InterSpec.

- A GEANT4 installation is **48 MB of libraries plus 2.1 GB of data**
  (726 MB of it just the low-energy EM data set).
- General-purpose transport cannot assume anything about your problem, so it
  cannot apply the variance reduction that a detector-efficiency calculation
  invites — the source is small, the detector subtends a known solid angle, and
  the quantity of interest is a peak area.

CeeLo is the narrow-purpose alternative: **a 1.8 MB static library (Release)
with every cross-section table compiled in and no runtime data files at all**,
built around
a biasing policy written specifically for this problem. It is fast enough to run
interactively inside a spectroscopy tool, which is where it is used.

The trade is deliberate: CeeLo does photon transport in the 30 keV – 3 MeV
regime and nothing else. No neutrons, no hadronic physics, no optical photons,
no beamline.

---

## Relationship to InterSpec

CeeLo is developed as part of [InterSpec](https://github.com/sandialabs/InterSpec),
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

The test suite takes ~15–30 minutes: these are Monte Carlo simulations with
statistical convergence checks, not unit tests of arithmetic.

**Dependencies:** C++17 compiler, CMake ≥ 3.15, Eigen3, and Boost (Boost.Test,
for the tests only).

**Optional features**, all off by default so that a plain checkout builds with
nothing extra:

| CMake option | Default | What it adds |
|---|---|---|
| `CEELO_WITH_SANDIADECAY` | OFF | Cascade-summing adapter driven by real decay data ([SandiaDecay](https://github.com/sandialabs/SandiaDecay)) |
| `CEELO_BUILD_RESPONSE` | OFF | Storable, MC-parameterized detector-response payload + EFFTRAN-style efficiency transfer (needs header-only RapidXml, deliberately not vendored) |
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
    // 3"x3" NaI(Tl): radius 3.81 cm, length 7.62 cm.
    calc.set_detector(DetectorShape::Cylinder, &nai, {3.81, 7.62});
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

## How it was developed

CeeLo is validated against **GEANT4 11.4.0** (`FTFP_BERT` +
`G4EmStandardPhysics_option4`), and the comparison is set up so the two codes
cannot silently disagree about the problem: **the GDML geometry and the GPS
source macro that GEANT4 runs are exported by CeeLo itself**, so geometry and
source are identical by construction rather than by careful transcription.

The GEANT4 results are committed in
[`tests/data/geant4_reference/`](tests/data/geant4_reference/) alongside CeeLo's
own in [`tests/data/ceelo_reference/`](tests/data/ceelo_reference/), so anyone
can re-run the comparison without installing GEANT4:

```bash
python3 profiling/compare_validation.py
```

### Validated configurations

Agreement of CeeLo against 16M–64M-event GEANT4 references, 30 keV – 3 MeV:

| # | Detector | Attenuator | Source | FEP | Total |
|---|---|---|---|---|---|
| 1 | 3″×3″ NaI | bare | point, 10 cm on-axis | ≤ 0.9 % | ≤ 0.6 % |
| 2 | 3″×3″ NaI | 1 mm Al | point, 10 cm on-axis | ≤ 0.7 % (≤ 1173 keV); −2.1 % @ 3000 | ≤ 0.6 % |
| 3 | 2″×2″ LaBr₃ | 0.5 mm Al | point, 5 cm on-axis | ≤ 0.6 % | ≤ 0.5 % |
| 5 | 1×1×0.5 cm CZT | bare | point, 5 cm on-axis | ≤ 1.3 % (30–662 keV); −5.3 % @ 1500 | ≤ 0.5 % |
| 6 | 3″×3″ NaI | bare | point, 15 cm, 45° off-axis | ≤ 0.6 % | ≤ 0.6 % |
| 7 | 3″×3″ NaI | 1 mm Al + 2 mm Pb | point, 15 cm on-axis | ≤ 0.8 % (≥ 200 keV) | ≤ 0.9 % (≥ 200 keV) |
| 8 | 3″×3″ NaI | 0.5 mm Al | Marinelli beaker, water | ≤ 1.2 % (≥ 100 keV); −2.3 % @ 59 | ≤ 1.1 % |
| 11 | 3″×3″ NaI | bare | point, 10 cm + 0.5 cm Fe shield | ≤ 0.6 % (200–2000 keV); −1.7 % @ 3000 | ≤ 0.6 % |
| 12 | 3″×3″ NaI | bare | 10×15×20 cm SS304 box, cellulose | ≤ 0.2 % (≥ 200 keV); −3.8 % @ 59 | ≤ 1.4 % (≥ 200 keV); −3.0 % @ 59 |

Measured against the committed GEANT4 references at ~0.3 % Monte Carlo
precision, so a few tenths of a percent of every entry is statistical.

Two families account for every entry above ~1 %, and both are documented under
[Shortcomings](#shortcomings): thin-crystal electron escape (config 5 above
1 MeV) and the high-energy bremsstrahlung/CSDA approximation (configs 2 and 11
at 3000 keV). The 59 keV entries for configs 8 and 12 are low-energy
scattered-photon effects in extended, self-attenuating sources.

Note that CeeLo's photon data is **EPICS2023**, while GEANT4 11.4.0 ships
**EPICS2017**. Where the two evaluations differ, disagreement with GEANT4 is
expected and is not by itself evidence of an error in either code.

The clearest case is **config 12 at 59 keV** (−3.8 % FEP, −3.0 % total), which is
almost entirely a data difference: iron's photoelectric cross section differs
between the two evaluations by ~2 % around 58 keV, and the source box's 2 mm
steel wall has an optical depth of ~2 there, so ~2 % in μ becomes ~4 % in
transmitted flux. Checked against the NIST *X-Ray Mass Attenuation Coefficients*
tables — an independent standard — **CeeLo is the closer of the two**: at
50–60 keV EPICS2023 sits +0.3 % from NIST for Fe, Cr and Ni, while the
EPICS2017-era values are ~1.3 % low. The effect falls to ≤0.7 % by 100 keV and
vanishes by 200 keV.

One documented exclusion: **config 7 at 100 keV**, where the ~−25 % total
discrepancy was traced to a GEANT4-side artifact (its in-tracking μ_total for Pb
is smeared across the 88 keV K-edge relative to its own exported tables, which
match NIST XCOM to 0.23 %). It is skipped in the automated gate.

Spherical and annular source geometries, and the cascade-summing path, are
validated separately — see [DESIGN.md](DESIGN.md).

---

## Performance compared to GEANT4

Both codes were run on the **same exported geometry and source**, both
cone-biased, both stopping when the **FEP efficiency reached 0.3 % relative
statistical precision** — so this is time-to-equal-answer, not time-per-event.
Single-threaded, Apple M1 Pro, Apple clang 21, `-O2 -march=native`.

| Config | Energy | CeeLo CPU | GEANT4 CPU | Speed-up |
|---|---|---|---|---|
| 1 — NaI 3″×3″, bare | 59 keV | 0.05 s | 2.11 s | **43×** |
| | 662 keV | 0.92 s | 8.16 s | **8.9×** |
| | 1332 keV | 2.14 s | 15.30 s | **7.2×** |
| | 2614 keV | 5.77 s | 33.58 s | **5.8×** |
| 3 — LaBr₃ 2″×2″ + Al | 59 keV | 0.05 s | 2.16 s | **45×** |
| | 662 keV | 1.26 s | 10.41 s | **8.2×** |
| | 1332 keV | 2.90 s | 19.55 s | **6.7×** |
| | 2614 keV | 7.90 s | 43.17 s | **5.5×** |
| 5 — CZT 1×1×0.5 cm | 59 keV | 0.07 s | 3.48 s | **49×** |
| | 662 keV | 2.59 s | 23.18 s | **9.0×** |
| 7 — NaI + Al + 2 mm Pb | 200 keV | 1.25 s | 81.78 s | **66×** |
| | 662 keV | 2.46 s | 23.03 s | **9.4×** |
| | 1332 keV | 4.76 s | 37.61 s | **7.9×** |
| | 2614 keV | 12.61 s | 81.48 s | **6.5×** |

**Median speed-up 8.6×; 5.5–9.4× across most of the range, 43–66× at low energy
and behind shielding.** Subtracting GEANT4's ~1.0 s fixed initialization cost the
median is still 7.7×. (Measured Aug 2026 on the EPICS2023 data.)

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
| Library | `libCeeLo.a` — **1.8 MB** (`Release`) | `lib/` — **48 MB** |
| Runtime data files | **none** | `G4EMLOW` — **726 MB** (EM only; full data set 2.1 GB) |
| Physics data, compiled in | **~575 KiB** of constant data | n/a — loaded from disk |
| Self-contained example binary | **~970 KB** stripped | — |

Library size is quoted for a `Release` build. The default build type is
`RelWithDebInfo`, where the static archive carries debug symbols and is ~27 MB;
those symbols do not reach a linked, stripped binary.

CeeLo has no runtime data directory, no environment variables to set, and
nothing to install alongside the binary.

---

## Where the data comes from

All physics data is compiled into the library. Full citations, licences, and
regeneration instructions are in
[THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).

| Data | Source |
|---|---|
| Photon cross-sections (photoelectric, K-shell photoelectric, Compton, Rayleigh, pair) | Direct **EPICS2023 EPDL** records, with process-specific grids retaining thresholds and absorption-edge sides |
| Incoherent scattering functions S(x,Z), coherent Rayleigh sampling CDFs | Direct **EPICS2023 EPDL**; form factors are used during offline CDF generation but are not embedded as runtime lookup tables; Compton-profile support via BSD-licensed [xraylib](https://github.com/tschoonj/xraylib) 4.2.1 |
| Bremsstrahlung spectra (Seltzer–Berger) | Public-domain **NIST EPQ** at pinned commit `ce1e589`, with one fixed low-κ interpolation rule |
| Atomic relaxation: K and L fluorescence yields, lines, Coster–Kronig (Z=1–99 for decay daughters) | Direct **EPICS2023 EADL** MF=28/MT=533 |
| Electron stopping powers | Direct **NIST ESTAR** collision and radiative values for every Z=1–92 at all 53 nodes from 10 keV–20 MeV |
| Material compositions | **ICRU Report 37** (air), **IAEA-TECDOC-1011**, CRC Handbook |
| Nuclear decay data (cascade summing, optional) | **ENSDF** via [SandiaDecay](https://github.com/sandialabs/SandiaDecay) |

The bremsstrahlung generator constructs the Geant4-aligned 32-node κ grid from
NIST EPQ uniformly: κ=0.025 uses a quadratic through {0, 0.05, 0.1}, κ=0.075
uses a cubic through {0, 0.05, 0.1, 0.15}, and all remaining nodes are copied
from equal-κ EPQ values. This is general interpolation of public-domain data,
not extraction or fitted reproduction of G4EMLOW coefficients. Aggregate
Geant4 spectra/rates remain committed as validation observables.

---

## Shortcomings

**Scope**
- Photon transport only, roughly 30 keV – 3 MeV. No photonuclear reactions, no
  neutrons, no optical photon transport, no hadronic physics.
- Only the geometries in the validated-configuration table have been checked
  against GEANT4. Everything else is extrapolation.
- Air between the source and the detector is **not transported**. A deterministic
  correction is planned; for now, long air paths are unmodelled.

**Known accuracy gaps**
- **Source-material fluorescence is not emitted.** After photoelectric
  absorption in the *source* material or its shielding, no characteristic X-ray
  is produced. **Total** efficiency is therefore low for high-Z self-attenuating
  sources — a thorium sphere runs −22 % / −4 % / −2.6 % in total efficiency at
  238.6 / 583 / 911 keV against GEANT4. FEP is unaffected (the X-ray cannot land
  in the primary peak), and low-Z sources are unaffected.
- **High-energy FEP** agrees to ≈0.2 % pooled above 2 MeV (11 gate rows at
  ~0.1 % MC precision) after the Aug 2026 electron-transport fixes
  (path-consistent Highland scattering + per-step energy-loss straggling in
  the crystal walk). The largest remaining single-config residual is LaBr₃ at
  +0.5–0.6 % ≥ 1 MeV.
- **Thin-crystal electron escape**: CZT FEP, formerly 2–5 % low above
  800 keV, now agrees within statistics (+0.7 / −0.3 / −1.3 % at
  800 / 1000 / 1500 keV, |z| < 0.7) after the same fixes.
- **Low-energy extended sources**: a water Marinelli beaker runs −2.3 % FEP at
  59 keV — small-angle Compton and forward Rayleigh that should land inside the
  peak window are under-recovered. It shrinks to ≤1.2 % by 100 keV.
- **Source-electron skin escape** uses an empirical, GEANT4-calibrated
  exit-state albedo gate rather than a full Mott / Goudsmit-Saunderson multiple
  scattering treatment. It is validated across Z = 5–82 to 1–2 %, but it is a
  calibration, not first principles. Two first-principles alternatives are
  selectable at compile time (`CEELO_SOURCE_ESCAPE_MODEL=tails|gs`) and are
  currently slightly less accurate on box geometries.
- **Efficiency transfer** (`CEELO_BUILD_RESPONSE`) is ≤ 1 % on-axis in the far
  field but a few percent near-field and off-axis; queries carry inflated
  uncertainties and provenance flags in those regimes.

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
  io/             Storable detector response, efficiency transfer (optional)
  materials/      Material definitions and built-ins
  physics/        Electron CSDA, stopping powers, multiple scattering
  transport/      Photon transport: Compton, PE, Rayleigh, pair, fluorescence
examples/       Runnable examples and benchmark configurations
tests/          Boost.Test suites (~35 ctest targets) + reference data
tools/          Cross-section regeneration, GEANT4 validation harness
profiling/      The GEANT4 validation gate script
visualization/  Dependency-free WebGL viewer for exported GDML geometries
```

Further reading:

- **[DESIGN.md](DESIGN.md)** — physics choices, biasing policy, geometry model,
  the full GEANT4 validation procedure, and detailed known limitations.
- **[TODO.md](TODO.md)** — the open-items register.
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
