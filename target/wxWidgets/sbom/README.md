# InterSpec Windows (wxWidgets) SBOM

This directory generates a [CycloneDX](https://cyclonedx.org/) Software Bill of Materials
for the **Windows wxWidgets desktop build** of InterSpec. The SBOM (`InterSpec_sbom.cdx.json`)
is produced by CI, written into the install directory, and therefore ships **inside the
Windows `.zip`** attached to the "Development Build" GitHub release.

## Why a curated list instead of a scanner

InterSpec has **no package manifest** (no vcpkg/Conan/npm-for-C++). Its Windows dependencies
come from three places:

1. Prebuilt static libraries built by `target/patches/dep_build_msvc2022.bat` (pinned versions).
2. Git submodules compiled in-tree (SpecUtils, SandiaDecay, QR-Code-generator).
3. Vendored source under `external_libs/` (Minuit2, Cuba, muParserX, pugixml, Faddeeva, and
   SpecUtils' header-only `date` / `nlohmann-json` / `rapidxml`).

Every major dependency is **statically linked** into `InterSpec.exe` / `InterSpec.dll`. So:

- A **binary scanner** (e.g. syft) finds nothing — there are no dynamic deps to enumerate.
- A **source-tree scanner** over-reports badly, pulling in things that are *not* in the
  Windows build (see the exclusion list below).

The accurate, low-noise result is a hand-curated component list (`components.json`) that a
small generator (`make_sbom.py`) turns into standard CycloneDX, plus a CI drift-guard.

## Files

- **`components.json`** — the single source of truth. Human-editable list of exactly the
  components in the final Windows build. Each entry has a `version_source` tag that ties it to
  an authoritative source and drives the drift check.
- **`make_sbom.py`** — pure-stdlib Python 3 generator. `generate` (default) writes the SBOM;
  `--check` validates `components.json` against the sources.

## Usage

```sh
# Generate the SBOM
python target/wxWidgets/sbom/make_sbom.py \
  --repo . \
  --git-sha "$(git rev-parse HEAD)" \
  --timestamp "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  --output InterSpec_sbom.cdx.json

# Drift guard (run in CI; fails the build if the list is stale)
python target/wxWidgets/sbom/make_sbom.py --check --repo .
```

The generator resolves submodule commit SHAs from the checked-out tree, so run it **after**
`git submodule update --init` (CI checks out with `submodules: true`).

## What is INCLUDED

The InterSpec application itself (top-level component, with version + build git SHA), plus:
Boost, zlib, Wt, Eigen, Ceres Solver, wxWidgets, Microsoft WebView2 Loader; SpecUtils,
SandiaDecay, QR-Code-generator (C++ variant); Minuit2, Cuba, muParserX, pugixml, Faddeeva;
the header-only `date`, `nlohmann-json`, `rapidxml`; and the D3.js v3 runtime asset shipped
in the installer.

## What is deliberately EXCLUDED (and why)

These live in the source tree but are **not** part of the Windows build, so they are absent
from the SBOM by design:

- **QR-Code-generator non-C++ variants** — `c/`, `java/`, `java-fast/`, `python/`, `rust/`,
  `rust-no-heap/`, `typescript-javascript/`: only `cpp/qrcodegen.cpp` is compiled.
- **Other platform targets** — `target/ios`, `target/android`, `target/electron`,
  `target/osx`, `target/macOsQuickLook`, `target/docker`, `target/sharedlib` (and the
  `ios-cmake` submodule): not used by the wxWidgets Windows app.
- **Dev/test tooling** — `target/testing`, `target/fuzzing`, `target/peak_fit_improve`,
  `target/example_code`, `target/idb_enrichment_check`; per-library `test`/`examples`/`tools`
  dirs; the `.venv_fission` Python virtualenv (incl. its PowerShell `Activate.ps1`).
- **SpecUtils Python bindings** (`SpecUtils_PYTHON_BINDINGS=OFF`) and **inja** templates
  (`SpecUtils_INJA_TEMPLATES=OFF`).
- **Compiled-out features** — MySQL, OpenSSL/TLS, FastCGI, Postgres (Wt built with these OFF),
  and libharu/PDF/libpng (Apple-only path).

## Updating after a dependency bump

When a dependency version changes (edit `dep_build_msvc2022.bat`, move a submodule, or update
a vendored dir), update the matching entry in `components.json`. The CI **drift guard**
(`--check`) fails the build if `components.json` no longer matches the pinned tokens in
`dep_build_msvc2022.bat` or if a listed submodule/vendored path is missing — so a stale SBOM
cannot ship silently.

Note: license identifiers and the versions of vendored/header-only libraries are best-effort
SPDX values; confirm against each upstream `LICENSE` when bumping.
