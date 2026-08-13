# Cross-section preparation

This directory regenerates CeeLo's compiled physics tables from locked,
permissively distributed primary sources. Normal CMake builds are offline: the
upstream archives are used only in a gitignored preparation cache and are never
vendored into the library.

The authoritative input record is [`sources.lock.json`](sources.lock.json). It
contains source URLs, versions, hashes, access dates, terms, citations, and every
production transformation. The retained notices are in [`source_notices/`](source_notices/).

## Production sources

| Runtime data | Locked source | Generator |
|---|---|---|
| Bremsstrahlung spectra | NIST EPQ commit `ce1e5895892b17854d857e596250e4ae70d6b274` | `sb_tables_epq.py`, `generate_element_support.py` |
| Photon cross-sections and angular factors | direct EPICS2023 EPDL | `generate_epics_photon.py` |
| K/L atomic relaxation | direct EPICS2023 EADL, Z=1–99 for decay daughters | `generate_eadl_relaxation.py` |
| Collision and radiative stopping | NIST ESTAR, all Z=1–92 | `generate_estar_stopping.py` |
| Compton-profile support and atomic weights | xraylib 4.2.1 (BSD-3-Clause) | `generate_element_support.py` |

No production generator reads Geant4 or G4EMLOW. Geant4 is used only by the
separate observable-validation tooling under `tools/geant4_validation/`.

## Prepare and verify inputs

From the repository root:

```bash
python3 tools/prepare_cross_sections/fetch_sources.py
```

This downloads into `tools/prepare_cross_sections/data/source-cache/`, verifies
each archive or logical table digest, and safely extracts archives. Re-running
with `--offline` proves that a prepared cache is complete:

```bash
python3 tools/prepare_cross_sections/fetch_sources.py --offline
```

The xraylib generator must run with xraylib 4.2.1 importable. The source archive
itself is also locked and verified so the BSD notice and exact source version are
recorded independently of the local Python installation.

## Regenerate

Set paths to the verified cache outputs (the fetch command prints each one):

```bash
CACHE=tools/prepare_cross_sections/data/source-cache
EPQ="$CACHE/extracted/nist_epq/EPQ-ce1e5895892b17854d857e596250e4ae70d6b274/src/gov/nist/microanalysis/EPQLibrary/BergerSeltzerBrem"

python3 tools/prepare_cross_sections/generate_element_support.py \
  --epq-tables "$EPQ" \
  --xraylib-archive "$CACHE/downloads/xraylib-4.2.1.tar.gz"
python3 tools/prepare_cross_sections/sb_tables_epq.py \
  --epq "$EPQ" --output src/cross_sections/element_data.cpp
python3 tools/prepare_cross_sections/generate_epics_photon.py \
  --epdl "$CACHE/downloads/EPDL2023.ALL"
python3 tools/prepare_cross_sections/generate_eadl_relaxation.py \
  --eadl "$CACHE/downloads/EADL2023.ALL"
python3 tools/prepare_cross_sections/generate_estar_stopping.py \
  --estar-dir "$CACHE/generated/nist_estar_text"
python3 tools/geant4_validation/validate_photon_migration.py
python3 tools/geant4_validation/validate_stopping_migration.py \
  --estar-dir "$CACHE/generated/nist_estar_text"
cmake --build build-rel -j8 --target CeeLo
python3 tools/prepare_cross_sections/coefficient_accounting.py --build-dir build-rel
python3 tools/prepare_cross_sections/update_generated_manifest.py
python3 tools/prepare_cross_sections/check_provenance.py
```

Run `update_generated_manifest.py` last. CI recomputes its hashes and fails if a
generator, source lock, generated table, or coefficient/error report changes
without a complete regeneration.

## Bremsstrahlung transformation

EPQ supplies the same 57 native incident-energy nodes used by the related
Seltzer–Berger family. CeeLo retains its 27 nodes from 10 keV through 20 MeV and
the current 32-node reduced-photon-energy grid. For every Z and every energy:

- κ=0.025 is quadratic interpolation through EPQ κ={0, 0.05, 0.1}.
- κ=0.075 is cubic interpolation through EPQ κ={0, 0.05, 0.1, 0.15}.
- Every other output node is copied from the equal-κ EPQ node; EPQ κ=0.125 is omitted.
- Each element is packed to `uint16` with one float32 decode scale. No row
  normalization is applied; the generation report gates packing error in row
  shape, integrated rate, and mean photon fraction.

EPQ `phi_rad` is an independent first-moment checksum, not a normalization
target. The transformation is general interpolation of public-domain EPQ data;
it has no element- or energy-specific correction and is not a transcription or
fit of G4EMLOW coefficients.

## Accuracy and size gates

Photon curves have process-specific descriptors over the historical 10 keV–10
MeV range. Endpoints, thresholds, and both sides of discontinuities are always
retained. Every other native EPDL node remains a validation probe, but is stored
only when adaptive refinement needs it; 16 logarithmic interior probes are also
tested in every interval. Values use per-curve uint16 packing. Rayleigh total
cross sections jointly refine one float32 energy grid for each four adjacent
elements, with every element independently required to pass the gates. Other
photon-process energy coordinates use exact uint16 indices into a fully
deduplicated shared float32 pool, so absorption-edge flanks stay distinct;
compact descriptors address flat streams by uint16 offset instead of retaining
per-curve pointers. Raw decoded photon-grid views are not retained. Refinement
continues until maximum material error is at most 0.1% (0.05% for pair
production) and integrated error is at most 0.05%. Runtime incoherent
scattering-function data starts at
0.01 inverse angstrom. Form-factor curves are reconstructed and validated only
inside the generator; runtime coherent scattering uses their offline CDF and no
form-factor lookup table is emitted. Angular behavior is gated against the
underlying lower-x EPDL data at 0.25% CDF error and 2×10⁻³ mean-cosine error.
K/L relaxation lines use compact variable-length arrays.
Each 53-node ESTAR stopping row uses a per-element affine `uint16` log-domain
encoding with float32 offset and scale; its native-node decode error is reported
separately from the upstream table precision.
The 256-node Rayleigh inverse CDF is generated offline and stored as a normalized
`uint16` row per element on one shared float32 grid. Quantization is gated at
8×10⁻⁶ absolute CDF error and the complete representation is gated against a
dense reference integration at 0.25% CDF error. No runtime CDF construction or
per-element CDF allocation is required.
Accuracy takes precedence over memory targets.

Generated CSV, JSON, and Markdown reports live in [`reports/`](reports/).
`coefficient_sizes.*` records separate baselines for `main`,
`feature/FullSpecApprox`, and `feature/env-scatter`, including native/retained
counts, encodings, grids, metadata, compiled constant bytes, and interpolation
errors. The three before-states are pinned by commit hash in the accounting
generator; shared grids are reported separately.

Numerical migration reports preserve measured failures rather than rewriting
them after review. A maintainer may accept an understood source-representation
or historical-approximation change in `reports/migration_acceptance.json`. Each
decision is bound to the exact SHA-256 of the failed report, so regeneration or
any changed result makes the acceptance stale and fails the provenance gate.

## Manual bremsstrahlung merge gate

The aggregate-only review generator may read a developer's local G4 data and
historical git objects, but writes no raw G4 cross-section rows:

```bash
python3 tools/geant4_validation/generate_brems_review_bundle.py \
  --g4-brem-dir "$G4LEDATA/brem_SB" \
  --estar-dir "$CACHE/generated/nist_estar_text"
```

Open `validation/bremsstrahlung_migration_review/review.html`, inspect the
specified element/material/energy cases, then record the decision in
`manual_review.json`. To enforce that gate locally or in a branch migration:

```bash
python3 tools/prepare_cross_sections/check_provenance.py --require-manual-approval
```

Aggregate Geant4 spectra, rates, efficiencies, and histograms are valid
validation fixtures when their configuration and statistical provenance are
recorded. Raw G4 cross-section columns must not be committed.

## Tests

```bash
python3 tools/prepare_cross_sections/test_sb_tables_epq.py
cmake --build build-rel -j8
ctest --test-dir build-rel --output-on-failure
```
