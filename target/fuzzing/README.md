# InterSpec libFuzzer harnesses

A small set of libFuzzer-based fuzz targets that exercise the parsers most exposed
to untrusted input: app URLs, JSON color themes, GADRAS DRF CSV, RelAct CSV, and
state XML.  Built with `-fsanitize=fuzzer,address,undefined` so a crash, ASan
report, or UBSan report fails the run.

This is a **standalone CMake project** — separate from the unit-test build.  The
sanitizer flags applied to `InterSpecLib` here would prevent the test binaries
from linking, so the two builds are kept in different directories.  Run these
periodically before tagging a release, or after touching any of the parsers
above.

## Toolchain — must use Homebrew LLVM

The Apple-distributed Xcode `clang` does not ship the libFuzzer runtime needed
for `-fsanitize=fuzzer` combined with ASan/UBSan.  Use the Homebrew LLVM
toolchain (the same one InterSpec's `external_libs/SpecUtils` fuzzers use):

```bash
brew install llvm
```

## Configure & build

From the repository root:

```bash
mkdir -p target/fuzzing/build
cd target/fuzzing/build

cmake -DCMAKE_C_COMPILER=$(brew --prefix llvm)/bin/clang \
      -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++ \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_PREFIX_PATH=/Users/wcjohns/install/wt-3.7.1_cxx20_prefix \
      -DWt_INCLUDE_DIR=/Users/wcjohns/install/wt-3.7.1_cxx20_prefix/include \
      -DBoost_DIR=/Users/wcjohns/install/wt-3.7.1_cxx20_prefix/lib/cmake/Boost-1.84.0 \
      -DEigen3_DIR=/Users/wcjohns/install/wt-3.7.1_cxx20_prefix/share/eigen3/cmake \
      -DCeres_DIR=/Users/wcjohns/install/wt-3.7.1_cxx20_prefix/lib/cmake/Ceres \
      -G Ninja \
      ..

cmake --build .
```

Adjust the `-D...DIR=` paths for your Wt prefix.  The defaults above match the
maintainer's local layout.

## Run

Use the `run_fuzzer.sh` wrapper at the project root.  It copies the seed corpus
to a temp directory so libFuzzer's coverage-expansion writes don't pollute the
committed seeds, sets sensible `ASAN_OPTIONS` defaults, and saves any crash
artifacts to `target/fuzzing/crashes/<fuzzer-name>/` before cleaning up.

```bash
cd target/fuzzing
./run_fuzzer.sh fuzz_color_theme -max_total_time=300
./run_fuzzer.sh fuzz_app_url -max_total_time=600 -rss_limit_mb=4096
```

The wrapper also passes `-print_final_stats=1` so each run reports coverage
and corpus growth.

If you want to run the binaries directly without the wrapper, do this so the
seed corpus stays clean:

```bash
WORK=$(mktemp -d)
ASAN_OPTIONS=detect_container_overflow=0:detect_leaks=0:abort_on_error=1 \
  ./build/fuzz_color_theme -max_total_time=300 -print_final_stats=1 \
                           -artifact_prefix="${WORK}/" \
                           "${WORK}" ../fuzzing/corpus/color_theme/
# move any crash files out of $WORK before:
rm -rf "${WORK}"
```

Why two corpus dirs: libFuzzer reads seeds from every directory listed but only
*writes* new coverage-expanding inputs to the **first** directory.  Pointing at
a temp dir first keeps the committed seed corpus untouched.

The recommended `ASAN_OPTIONS`:

- `detect_container_overflow=0` silences a known libc++/libFuzzer false-positive
  in libFuzzer's own (un-instrumented) corpus loader.
- `detect_leaks=0` ignores per-iteration leaks (`fuzz_app_url` leaks every Wt
  widget by design across iterations; `-rss_limit_mb` is the right knob there).
- `abort_on_error=1` makes libFuzzer record a crash file when ASan/UBSan fires.

Useful flags:

- `-max_total_time=N` — run for N seconds
- `-runs=N` — run N iterations then exit (`-runs=0` = single replay of corpus)
- `-print_final_stats=1` — print coverage / corpus size stats on exit
- `-jobs=4 -workers=4` — parallel fuzzing
- `-rss_limit_mb=4096` — kill if memory grows past 4 GB (recommended for fuzz_app_url)
- `-only_ascii=1` — restrict mutations to printable ASCII (useful for JSON / CSV / URL fuzzers)

If a crash is found, libFuzzer writes the offending input to a file named
`crash-<sha1>` in the working directory.  Replay it with:

```bash
./fuzz_color_theme crash-deadbeef...
```

## Targets

- `fuzz_app_url` — full `interspec://...` URL dispatch via
  `InterSpec::handleAppUrl`.  Builds a headless `InterSpec` instance once at
  startup using `Wt::Test::WTestEnvironment` and the static-data dir located
  via `INTERSPEC_DATA_DIR` or by walking up from the cwd.  Each fuzz input
  becomes one URL, so this target exercises every per-tool `handleAppUrl`
  (decay, simple-mda, dose, flux, gammaxs, drf, etc.) through the same path
  the `EnterAppUrlWindow` uses at runtime.  Memory-state accumulates across
  iterations - run with `-rss_limit_mb=2048` (or similar) to bound it.
- `fuzz_color_theme` — `ColorTheme::fromJson`
- `fuzz_rel_eff_drf_csv` — `DetectorPeakResponse::parseSingleCsvLineRelEffDrf`,
  the single-line relative-efficiency DRF CSV format used in `data/common_drfs.tsv`.
- `fuzz_gadras_efficiency_csv` — `DetectorPeakResponse::fromGadrasDefinition`
  with a fixed known-good `Detector.dat`; mutates the `Efficiency.csv` side.
  Seed corpus is `data/GenericGadrasDetectors/<detector>/Efficiency.csv`.
- `fuzz_gadras_detector_dat` — `DetectorPeakResponse::fromGadrasDefinition`
  with a fixed known-good `Efficiency.csv`; mutates the `Detector.dat` side.
  Exercises both the old-style line-numbered parser and the XML branch
  (selected when the first line contains "xml").  Seed corpus is
  `data/GenericGadrasDetectors/<detector>/Detector.dat` plus a small synthetic
  XML example.
- `fuzz_rel_act_csv` — `RelActCalcManual::peak_csv_to_peaks`
- `fuzz_state_xml` — `SpecMeas` N42 file loader

## Suggested cadence

- After touching any parser: run the relevant fuzzer for ~5 minutes.
- Before a release: run each fuzzer for 30-60 minutes.
- After a major parser refactor: 8h+ overnight runs are worthwhile.


## Untargeted parsers — candidates for future fuzz harnesses

The list below is what an exhaustive grep across `InterSpec/` and `src/` turned
up: functions that read user-supplied bytes (files, URLs, JSON, XML) and that
are **not** currently fuzzed.  External libraries (`SpecUtils`, `SandiaDecay`)
are excluded — they have their own fuzzers.

Add a new harness whenever one of these functions sees significant changes, or
before any release if the function has changed since the last release.

### `interspec://` URL handlers — covered transitively, dedicated harnesses optional

`fuzz_app_url` dispatches through `InterSpec::handleAppUrl`, which transitively
exercises every per-tool handler (`DecayWindow::handleAppUrl`,
`DetectionLimitSimple::handleAppUrl`, `DoseCalcWidget::handleAppUrl`,
`ExportSpecFile::handleAppUrl`, `FluxTool::handleAppUrl`,
`GammaCountDialog::handleAppUrl`, `GammaXsGui::handleAppUrl`,
`OneOverR2Calc::handleAppUrl`, `RemoteRid::handleAppUrl`,
`SimpleActivityCalc::handleAppUrl`, `LicenseAndDisclaimersWindow::handleAppUrlPath`,
`DrfSelect::handle_app_url_drf`, etc.).  Adding a per-tool harness is only
worthwhile if a particular tool sees heavy churn or you want faster
iteration speed for that tool's specific corpus.

The DRF URL-encoded path used inside `parseSingleCsvLineRelEffDrf` is also
covered by `fuzz_rel_eff_drf_csv` (when fuzzed lines start with the
`UrlEncoded` field-3 marker).  `DetectorPeakResponse::fromAppUrl` (instance
member) and the `RADDATA://` / `interspec://G0/` spectrum URL paths are
**not** yet fuzzed and could use dedicated harnesses.

### Other DRF parsers (in `DetectorPeakResponse.{h,cpp}`)

Currently fuzzed: `parseSingleCsvLineRelEffDrf`, `fromGadrasDefinition`.  Not
fuzzed:

- `DetectorPeakResponse::parseMultipleRelEffDrfCsv` — multi-line variant of
  the rel-eff DRF CSV parser
- `DetectorPeakResponse::parseGammaQuantRelEffDrfCsv` — GammaQuant-format
  rel-eff CSV
- `DetectorPeakResponse::parseAngleOutxFile` — Angle Output X file
- `DetectorPeakResponse::parseEfficiencyCsvFile` — alternate Efficiency.csv
  parser used outside `fromGadrasDefinition`
- `DetectorPeakResponse::parseDetectorDatGeometry` — geometry-only parse
- `DetectorPeakResponse::fromXml` — DRF state restored from saved XML
- `DrfSelect::parseInterSpecRelEffCsvFile` — file-path entry point for the
  InterSpec rel-eff CSV format

### XML state deserializers

InterSpec saves analysis state into XML inside the N42 file (and as standalone
state files).  Loading a corrupted or malicious state file walks each of these
parsers.  None are currently fuzzed.

- `PeakDef::fromXml` and `PeakContinuum::fromXml` (in `PeakDef.h`)
- `PeakModel::*` peak-list deserialization paths
- `ReferenceLineInfo::deSerialize`
- `ReferencePhotopeakDisplay::deSerialize`
- `IsotopeSearchByEnergy::deSerialize` (string and node overloads), plus the
  `SearchEnergy::deSerialize` and category-list overloads
- `GammaInteractionCalc::ShieldingSourceFitOptions::deSerialize`
- `ShieldingSourceDisplay::deSerialize` (multiple overloads)
- `RelActAutoGui::deSerialize`
- `RelActManualGui::deSerialize`
- `RelActCalcAuto::*::fromXml` (RoiRange, NucInputInfo, FloatingPeak,
  RelEffCurveInput, Options, ManualActRatioConstraint, MassFractionConstraint,
  RelActAutoGuiState — about a dozen `fromXml` methods)
- `RelActCalcManual::*::fromXml` (ManualActRatioConstraint, etc.)
- `RelEffShieldWidget::fromXml`
- `SpecMeas::decodeSpecMeasStuffFromXml` (the InterSpec-specific N42 extension
  parsing — distinct from the SpecUtils-side N42 loader)

A single state-file fuzzer that feeds bytes into the top-level
`SpecMeas::load_from_N42` route would exercise most of these transitively;
direct per-class harnesses are useful when isolating regressions.

### JSON parsers (other than `ColorTheme::fromJson`)

- `EventXmlFilterInfo::parseJsonString` (in `SpecFileQueryDbCache.h`) — JSON
  used by the spectrum-file query tool; takes user-supplied filter definitions
- `FileDragUploadResource::handleRequest` — parses a JSON body containing a
  `fullpath` field (Electron / native-app drop path); already validated
  against path traversal but the parse step itself is unfuzzed
- `Wt::Json::parse` of the saved query in `SpecFileQueryWidget.cpp`
- `Wt::Json::parse` of `userOptionsJson` in `InterSpec.cpp` (saved user prefs)
- `Wt::Json::parse` of upload manifests in `InterSpecServer.cpp`
- `Json::parse` of the `controlInfo` JSON in `UseInfoWindow.cpp`
- `Json::parse` of remote-RID JSON responses in `RemoteRid.cpp` (input is
  trusted-ish since it comes from an executable the user configured, but
  worth fuzzing since the server can be remote)

### Spectrum file loaders (beyond N42)

`SpecMeas::load_from_N42` is fuzzed (see `fuzz_state_xml`).  Not fuzzed:

- `SpecMeas::load_from_iaea` — IAEA-format text spectrum (the InterSpec
  override that adds peak/state restoration on top of the SpecUtils base)
- `QRSpectrum::decode_spectrum_urls` and related — decoding spectra encoded
  into `RADDATA://` URLs / QR codes (InterSpec consumes these via
  `SpecMeasManager::handleNonSpectrumFile` and the multi-URL dialog)

### Peak CSV parsers

- `PeakModel::gadras_peak_csv_to_peaks` — InterSpec's own peak-CSV importer
  (distinct from `RelActCalcManual::peak_csv_to_peaks`, which **is** fuzzed)

### Internal data-file XML parsers

These read files shipped in `data/` so the threat model is mainly "did we
break something during refactoring" rather than active attack, but they're
worth coverage for crash regressions:

- `MoreNuclideInfo` XML loader (`data/sandia.MoreNuclideInfo.xml`)
- `ReactionGamma` constructor parsing `data/sandia.reactiongamma.xml`
- `XRayWidthServer` parsing `data/xray_widths.xml`
- `IsotopeSearchByEnergy` category-list XML

### `Wt::WResource` request handlers

Each `Wt::WResource` subclass exposes a URL endpoint; query parameters arrive
from arbitrary HTTP clients in web-deployment mode.  Only resources that
parse non-trivial input from the request need fuzzing:

- `FileDragUploadResource::handleRequest` (already noted under JSON above)
- `RestRidInputResource::handleRequest` (in `RemoteRid.cpp`)
- The various download resources (`PeakCsvResource`, `DecayCsvResource`,
  `DrfDownloadResource`, `JsonDownloadResource`, `DownloadSpectrumResource`,
  `DecayChainHtmlResource`) — these are output-only and don't parse input,
  so fuzzing has low value.

### How to add a new harness

1. Drop a new `fuzz_<thing>.cpp` in this directory with a single
   `extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)`.
2. Add `add_fuzz_target( fuzz_<thing> )` to `CMakeLists.txt`.
3. Create `corpus/<thing>/` with 3-10 representative inputs.
4. Add a one-line description to the "Targets" section above.

Use `gadras_baseline.h` as a template when your harness needs a fixed
known-good companion file (e.g., when fuzzing one half of a two-file format).
