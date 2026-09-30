# InterSpec Light

A browser-only, serverless subset of InterSpec, built as **one self-contained HTML file** that works
when opened straight from disk (`file://`), and makes no network requests (unless given a spectrum
URL).  Two variants are built:
- `dist/InterSpecLight.html` (~2.24 MB): the WASM module and reference lines are gzip'ed; needs a
  browser with `DecompressionStream` (Safari 16.4+, Chrome 80+, Firefox 113+).
- `dist/InterSpecLight_uncompressed.html` (~6.0 MB): plain base64 WASM and plain JSON, for older
  browsers, or email/web content scanners that object to the compressed data (it can resemble
  "HTML smuggling").

In both, the JavaScript is inline, minified but readable; nothing is decoded and then executed.

It can:
- show foreground, background, and secondary spectra, with background and secondary live-time normalized
- show a time chart for passthrough/search-mode files, where dragging selects the samples for each spectrum
- let you step through samples of multi-sample files, and show or hide detectors of multi-detector files
- show reference lines for ~290 nuclides at their default ages, NORM background, reactions, and element x-rays
- search one or two energies for nearby lines
- fit peaks: double-click, ctrl-drag, ROI-edge drag, shift-drag erase, and a right-click menu for delete, refit, continuum, and skew
- automatically search for and fit all peaks ("Search for peaks"; existing peaks are kept)
- assign sources to peaks: from the shown reference lines when fit, from the right-click menu's "Assign as"
  items, or by typing one into the peak table (e.g., "Cs137", "Pb 75", "Th232 S.E.")
- edit or fit (from peaks) the energy calibration, and apply it to the shown detectors
- export the energy calibration as a CALp file, or apply one (drop it on the page), as InterSpec does
- export N42/PCF/CHN/SPE/CNF/... files, and the peaks as an InterSpec-compatible peak CSV
- read and write peaks in N42-2012 and IAEA SPE files, the way InterSpec does (optional; see below)

The C++ is a hard fork of InterSpec's peak fitting and energy calibration code, plus SpecUtils, compiled
to WebAssembly with Emscripten; it uses Eigen and Ceres, but not Boost.  It is built separately from
InterSpec, and nothing in InterSpec depends on it.

## Building

### Prerequisites
- Git, CMake 3.20 or newer, Ninja, Python 3, and a C++20 compiler (for the native tools).  With
  Python's Pillow installed, the logos are embedded at their display size (otherwise at full size).
- [emsdk](https://emscripten.org/docs/getting_started/downloads.html), activated in the shell you build
  from; Emscripten 5.0.6 is what the page is built and tested with:
  ```bash
  git clone https://github.com/emscripten-core/emsdk.git
  cd emsdk && ./emsdk install 5.0.6 && ./emsdk activate 5.0.6
  source ./emsdk_env.sh
  ```
- Network access for the first build, which downloads the WASM build's dependencies at pinned versions
  (see `deps/build_deps.sh`): Eigen 5.0.1, Abseil LTS 2026-01, and Ceres (2026-03 master).  Boost is not
  needed.
- Optional: the dependency prefix InterSpec itself is built with (see `target/dep_build/`), given as
  `NATIVE_PREFIX`, to also build the native test driver `light_cli` (it needs Ceres).

### Build
From an activated emsdk shell:
```bash
cd target/InterSpec_light
./build.sh                                         # or: NATIVE_PREFIX=<prefix> ./build.sh
```
This writes `dist/InterSpecLight.html` and `dist/InterSpecLight_uncompressed.html`.  The first build
also downloads and builds the dependencies, which takes a few minutes (~2 minutes on a recent Mac with
a fast connection); the page itself then builds in about a minute.  All builds use 4 cores.
`dist/InterSpecLight.html` is tracked in git, so the page can be used without building it; a build
overwrites it, so commit it along with source changes.  The uncompressed variant is not tracked, to
avoid bloating the repository.

`./build.sh --without-peak-files` (CMake option `LIGHT_PEAK_FILE_IO=OFF`) leaves out reading and writing
peaks in N42 and SPE files, which saves ~27 kB of the compressed page (~70 kB uncompressed): the XML/CSV
code (`cpp/fork/src/PeakDefXml.cpp`, `cpp/app/PeakFileIo.cpp`), and each reference line's decay daughter
(which InterSpec needs to identify a peak's nuclear transition).

`build.sh` runs these steps (each can also be run on its own):
1. `deps/build_deps.sh` downloads Eigen, Abseil, and Ceres to `deps/src`, and cross-compiles
   them with Emscripten into `deps/prefix` (skipped once done).  To use other versions, change the
   commits in the script, and delete `deps/src`, `deps/build`, and `deps/prefix`.
2. `build_native/` builds `refgen` (and `light_cli`, if Ceres is found).
3. `refgen` writes `build_native/gen/ref_lines.json` from `data/sandia.decay.xml`, for the sources
   listed in `refgen/sources.txt` (seeded from `data/` by `refgen/collect_sources.py`; edit freely).
4. `build_wasm/` builds `light_wasm.js` + `light_wasm.wasm`, as MinSizeRel (`-Oz`), except the
   translation units with the hot fitting code (`PeakFitLM`, `PeakDists`, `PeakFit`), which stay at
   `-O3`: an all-`-Oz` build is ~3x slower to fit peaks.
5. `assemble_html.py` writes both pages.  The chart JS/CSS are read from InterSpec's sources at this
   step, and the JavaScript is minified with the terser that ships with Emscripten.  The WASM and
   reference lines are embedded as data blocks that `App.readAsset` (`web/main.js`) reads on load
   (~0.2 s either way).

## Using

Open `dist/InterSpecLight.html`, then drop a spectrum file on the page, or use an open button.  While
dragging a file, drop zones pick foreground, background, or secondary; each spectrum in the side panel
also has its own open button.  The Sandia logo at top right shows/hides the side panel.  The Help
section of the side panel lists the mouse gestures.

URL parameters: `?url=<spectrum>&background=<spectrum>&secondary=<spectrum>`.  When the page is
opened from `file://`, the server holding the spectra must send `Access-Control-Allow-Origin: *`.
Alternatively, serve the page and the spectra from the same web server.

The peak CSV is a port of InterSpec's `PeakModel::write_peak_csv` (the "Full" variant: same
headers, columns and formatting); InterSpec's own `PeakModel::csv_to_candidate_fit_peaks` reads it
back with sources, continuum/skew types, and colors intact.

### Peaks in N42 and SPE files (unless built `--without-peak-files`)
Files are read and written as InterSpec does, so either program opens the other's files with the
peaks, their sources, continua, and skews:
- **N42-2012** exports include InterSpec's `<DHS:InterSpec>` element: the peaks of every set of samples
  of the file, and which samples and detectors are displayed (InterSpec, and this page, show those
  samples on opening the file).  Files with this element keep their sample numbers on load, since the
  peaks are keyed by them (as InterSpec does).  Peaks InterSpec's automated search saved as fit hints
  are not loaded.
- **IAEA SPE** exports hold the displayed spectrum, with the displayed peaks in `$PEAK_INFO_CSV:` (the
  peak CSV) and `$PEAKLABELS:` sections.  Only the CSV is read back, so SPE peaks keep only what the
  CSV holds (values to its printed precision; no calibration or fit flags).
- Peak sources from SPE files are matched to this page's reference-line library; a line not in the
  library keeps its nuclide and energy, but not its decay, so InterSpec would not see it as assigned
  if the peak is later saved to N42.  Annihilation peaks are written without their positron decay
  (InterSpec accepts that).  As in InterSpec, values are read back as float (~7 digits).
- InterSpec itself mis-reads escape peaks from the peak CSV (and so SPE): it takes the CSV's escape-peak
  energy as the photon energy.  This page adds the 511/1022 keV back.

## Layout

| Path | Contents |
|---|---|
| `cpp/fork/` | **Hard fork** of InterSpec's peak fitting and energy-cal code: `PeakDef`, `PeakDists`, `PeakFit`, `PeakFitLM`, `PeakFitSpecImp`, `PeakFitUtils`, `PeakFitDetPrefs`, `EnergyCal`, and the `_imp.hpp` headers.  No Wt, SandiaDecay, Boost, or threads: `LightMath` has ports of the few Boost.Math functions used, and the Savitzky-Golay filter uses Eigen.  The peak XML (`PeakDefXml.cpp`) is only built with peak-file support.  Peak sources are a plain `PeakDef::Source` struct.  `DetectorPeakResponse.h` keeps only the FWHM forms; `LightWtShim.h` has minimal `WFlags`/`WColor`. |
| `cpp/app/` | `Session`: loaded files, display slots, sample/detector selection, peaks, energy calibration.  Ported from `SpecMeasManager`, `InterSpec`, `PeakSearchGuiUtils`, `D3SpectrumDisplayDiv`, `D3TimeChart`, and `EnergyCalTool`.  `RefLib`: the reference-line library, and peak sources (assigned when fit, suggested, or typed).  `PeakCsv`: the peak CSV writer.  `PeakFileIo`: peaks in N42 and SPE files (optional).  `LightApi.cpp`: the single `light_call(json)` entry point. |
| `cpp/tools/light_cli.cpp` | Native driver: runs JSONL request scripts through `light_call`, for tests, and debugging with lldb. |
| `deps/` | `build_deps.sh`, which downloads and builds the WASM build's dependencies. |
| `refgen/` | Reference-line generator, and its source list. |
| `web/` | Page template, CSS, and app JS.  All state lives in the WASM; each call returns the changed parts of the display. |
| `tests/` | `run_tests.sh` runs each `tests/*.jsonl` natively and through WASM (Node), compares the results, and checks the expectations scripts put on responses (`check_expectations.py`).  `browser_smoke.mjs` and `browser_interactions.mjs` drive the built page in headless Chrome with real mouse input.  `interspec_compat/` checks InterSpec itself reads the peak files this page writes; `data/` holds files InterSpec wrote (by `interspec_compat/run.sh --make-fixtures`). |

## Testing

```bash
NATIVE_PREFIX=<prefix> ./build.sh                   # light_cli is needed by run_tests.sh
tests/run_tests.sh                                  # native vs WASM parity, and expectations
node tests/browser_smoke.mjs <out_dir> [page]       # page, fitting, e-cal, export, search, time chart
node tests/browser_interactions.mjs <out_dir> [page] # drop zones, NaI, ctrl-drag, ROI-edge drag, erase
tests/interspec_compat/run.sh                       # InterSpec reads this page's N42/SPE peaks (optional)
```

Screenshots go to `<out_dir>`.  `[page]` defaults to `dist/InterSpecLight.html`; pass
`dist/InterSpecLight_uncompressed.html` to test that variant.  The browser tests use the Playwright
install in `target/testing/PlaywrightPhoneEmulation` (`npm install` there once), with Chrome from
`CHROME_PATH`, else its standard macOS location, else Playwright's own browser
(`npx playwright install chromium`).

`interspec_compat/run.sh` needs InterSpec itself built (`INTERSPEC_BUILD`: its build directory, default
`build_vscode`): it builds `interspec_check` against InterSpec's library, and checks InterSpec sees the
same peaks in N42 and SPE files this page wrote as in the ones InterSpec wrote.

## Simplifications relative to InterSpec
- There are no prompts on load.  Raw data is used when a file also has derived data, and the first
  energy binning when it has several.  Only "VD*" detectors are shown if a file has several of them.
- Peaks are fit only on the foreground.  There is no detector response, and no "hint" peaks from an
  automated search.
- The automated peak search is InterSpec's `ExperimentalAutomatedPeakSearch::search_for_peaks`, run
  serially with no progress or cancel (~0.5 s for a 16k-channel HPGe spectrum in Chrome; ~1 s for a
  busy Pu spectrum).  Its HPGe width filter uses a least-squares sqrt-polynomial FWHM fit to the found
  peaks in place of a detector response.  New peaks get sources from the shown reference lines,
  largest first.
- Energy-cal changes apply to all samples of the foreground file's shown detectors.  Other files
  are unchanged, and lower-channel-energy calibrations are read-only (though a CALp file can set or
  replace them).  Deviation pairs are kept as-is.
- A dropped CALp file is applied the same way (all samples, shown detectors), but with its own
  deviation pairs, which only reach the detectors that have the displayed calibration (others get the
  change propagated, as for an edit).  A calibration named for one detector is refused if the shown
  detectors have different calibrations, and a CALp file with per-detector calibrations needs one for
  each shown detector.  If a background or secondary from another file, with as many channels, is
  loaded, a dialog asks which spectra to apply it to.  "Revert" restores every loaded file's calibration.
- Peak sources come only from the reference-line library: typed sources are read as InterSpec's peak
  table reads them, but matched to the library's lines, so a nuclide not in the library can not be
  assigned.
- Reference lines use each source's default age only, with no shielding or detector effects.  The
  NORM set reuses InterSpec's precomputed soil-transport table.  To limit size, `refgen` keeps a line
  if its importance (intensity x sqrt(E)) is at least 0.1% of the source's largest, or if it is at
  least 10% of the largest importance within max(10 keV, 10% of E) (e.g., Pu239's 375 and 414 keV
  lines, which are ~1E-3 of its U L x-rays); lines below 3E-6 of the largest are always dropped.
- The WASM runs on the main thread (fits typically take well under a second).
