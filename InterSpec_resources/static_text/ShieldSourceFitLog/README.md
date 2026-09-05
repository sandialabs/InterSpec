# Activity/Shielding Fit & Batch Peak-Fit Report Template Authoring Guide

> **Intent.** This document is meant to be handed to an LLM together with a description of the
> report a user wants. The LLM should produce a complete Inja template (`*.tmplt.html`,
> `*.tmplt.txt`, `*.tmplt.csv`, `*.tmplt.md`, etc.) that consumes one of the JSON payloads
> documented below. Eight default templates ship in
> `InterSpec_resources/static_text/ShieldSourceFitLog/`; they are themselves examples of what
> the LLM should produce. The InterSpec code repository can be accessed at
> https://github.com/sandialabs/InterSpec/ , and Inja's documentation at
> https://github.com/pantor/inja / https://pantor.github.io/inja/ .
>
> **Two report families** are documented here:
> - **Activity / Shielding fit** — §5. The result of an Activity/Shielding source fit (one
>   spectrum at a time; or wrapped in a `Files[]` array for batch summaries). Use this family
>   for reports that show fitted activities, ages, shielding thicknesses, χ², per-peak
>   model-vs-observed comparisons, etc.
> - **Batch peak-fit** — §6. The result of an exemplar-driven peak fit applied to one or more
>   spectra (used by `--batch-peak-fit` and the "Batch Peak Fit" GUI). Reports list fitted
>   peaks, continua, optional energy-cal refit, and peak-source assignments.
>
> The two families share infrastructure (`foreground` / `background`, `EnergyCal`, chart
> assets, application metadata, the Inja env, the two custom callbacks `printFixed` /
> `printCompact`); see §5.4 / §5.5 / §5.14 / §5.15 / §7 / §8 / §9.
>
> **Last updated:** 2026-08-29

## 1. What this is

`InterSpec_resources/static_text/ShieldSourceFitLog/` is the directory in
[InterSpec](https://github.com/sandialabs/InterSpec) that holds [Inja](https://github.com/pantor/inja)
templates for two related tools:

- The **Activity/Shielding source fit** tool, which uses fitted gamma photo-peaks plus a
  user-supplied geometry (sphere / cylinder / rectangular, with optional shielding layers
  and self-attenuating or trace sources) to recover absolute nuclide activities, source
  ages, and shielding thicknesses. The C++ data model is `ShieldingSourceFitCalc::ModelFitResults`
  populated by `ShieldingSourceFitCalc::ShieldingSourceChi2Fcn` and friends. The JSON
  payload is built by `BatchInfoLog::shield_src_fit_results_to_json()` (see
  [`src/BatchInfoLog.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/BatchInfoLog.cpp)).
- The **Batch peak-fit** tool, which takes a user-curated "exemplar" peak list and refits
  those peaks in a batch of spectra. The JSON payload is built by
  `BatchInfoLog::add_peak_fit_results_to_json()` (same file) and per-file results are
  wrapped into `Files[]` for batch summaries.

Eight templates ship in `InterSpec_resources/static_text/ShieldSourceFitLog/` by default:

| File | Family | Purpose |
|---|---|---|
| `act_fit.tmplt.html`            | Activity / Shielding | Per-spectrum full HTML report with embedded D3 spectrum chart and shielding/source-fit plot |
| `std_fit_log.tmplt.txt`         | Activity / Shielding | Per-spectrum plain-text log |
| `std_summary.tmplt.html`        | Activity / Shielding | Multi-file batch-summary HTML (iterates `Files[]`) |
| `std_summary.tmplt.csv`         | Activity / Shielding | Multi-file batch-summary CSV |
| `peak_fit.tmplt.html`           | Batch peak-fit       | Per-spectrum full HTML peak-fit report |
| `std_peak_fit_log.tmplt.txt`    | Batch peak-fit       | Per-spectrum plain-text peak-fit log |
| `std_peak_fit_summary.tmplt.html` | Batch peak-fit     | Multi-file peak-fit summary HTML |
| `std_peak_fit_summary.tmplt.csv`  | Batch peak-fit     | Multi-file peak-fit summary CSV |

A user can drop their own `.tmplt.*` files into the **writable** template directory
(see §3) and InterSpec will pick them up automatically.

## 2. What these tools do

### 2.1 Activity / Shielding source fit

The tool takes a gamma-ray spectrum with peaks already fit and assigned to nuclides, plus
a user-supplied geometry, and recovers:

- **Nuclide activities** (`Sources[].Activity_*`) in absolute units (Bq / Ci variants).
- **Nuclide ages** when fittable (`Sources[].AgeIsFit`); InterSpec optionally couples ages
  for nuclides of the same element (`ActShieldFitSetup.FitOptions.ElementNuclidesSameAge`).
- **Shielding thicknesses and dimensions** for each shielding layer
  (`Shieldings.Shields[].Thicknesses_*`, `OuterDims_*`, `InnerDims_*`); shieldings can be
  either material-defined (chemical formula + density) or generic AN/AD.
- **Goodness of fit** (`FitChi2`, `NumDof`) and per-peak comparison
  (`PeakToModelComparison.UsedPeaks[]`).

Supported source kinds:

- **Point-like** (the default): activity sits at a single point at `Distance` from the detector.
- **Self-attenuating** (volumetric): the nuclide is distributed inside one of the shielding
  layers, attenuating its own emissions. Set via `Sources[].IsSelfAttenSource` /
  `Shieldings.Shields[].SelfAttenSources[]`.
- **Trace** (surface contamination, in-volume contamination, exponential depth distribution):
  `Sources[].IsTraceSource`, with `TraceActivityType` ∈ {`TotalActivity`, `ActivityPerCm3`,
  `ActivityPerGram`, `ExponentialDistribution`}.

Geometries: `Spherical`, `CylinderEndOn`, `CylinderSideOn`, `Rectangular`. The number of
shielding dimensions varies by geometry (1, 2, or 3); `Shieldings.NumDimensions` and
`Shieldings.DimensionMeanings[]` tell a template how to label them.

**Fixed-geometry detectors** (`FixedGeometryDetector = true`) hold a calibrated efficiency
appropriate for a specific source-detector configuration; in that mode there is no
`Distance` and `Geometry`, only `ActShieldFitSetup.FixedGeometryType` (a human-readable
description: `"total activity"`, `"activity per square centimeter"`, etc.).

### 2.2 Batch peak-fit

The batch peak-fit tool refits user-curated peaks in a series of spectra. The user
provides an **exemplar** (a previously-saved spectrum with peaks defined) and a list of
input files; for each input file the tool:

1. Loads the spectrum.
2. Optionally refits the energy calibration using a subset of exemplar peaks.
3. Refits each exemplar peak in the input spectrum.
4. Optionally subtracts a background spectrum's peaks.
5. Records which peaks did *not* converge (visible in templates as `NotFitPeaks.Peaks[]`)
   versus those that did (`FitPeaks.Peaks[]`).
6. Computes a detection limit ("MDA") for each peak that did not converge, at
   `NotFitPeaks.Peaks[].Mda` — see §6.6.1. Activity/shielding fits get the same section, with
   the limits also expressed as activities.

The exemplar's original peak set is also exposed for reference at `ExemplarPeaks.Peaks[]`.

## 3. Where templates live on disk

InterSpec looks for templates in two places: a **bundled** read-only directory (ships with
InterSpec) and a **per-user writable** directory (you control it).

### 3.1 Bundled (read-only) directory

```
<InterSpec install root>/InterSpec_resources/static_text/ShieldSourceFitLog/
```

This is where the eight default templates listed in §1 live. **Do not** edit files here:
they are overwritten on every InterSpec update, and on macOS / Windows the directory may
not even be writable. Copy any template you want to customise into the writable directory
(§3.2) first.

### 3.2 Per-user writable directory

InterSpec scans this directory at startup and dynamically adds any well-named templates
it finds to the "Calculation Log" dropdown in the Activity/Shielding GUI, and accepts
their filenames on the batch CLI:

```
<InterSpec::writableDataDirectory()>/act_shield_fit_templates/
```

The directory does **not** exist by default — you have to create it.

| OS | `<writableDataDirectory>` |
|---|---|
| **macOS**   | `~/Library/Application Support/sandia.InterSpec/` |
| **Windows** | `%APPDATA%\InterSpec\` (typically `C:\Users\<you>\AppData\Roaming\InterSpec\`) |
| **Linux (desktop / wxWidgets build)** | `$XDG_DATA_HOME/InterSpec/` if set, else `~/.local/share/InterSpec/` |
| **Linux (Electron build)** | `$XDG_CONFIG_HOME/InterSpecAddOn/` if set, else `~/.config/InterSpecAddOn/` |
| **Linux (batch CLI fallback)** | `~/.local/share/sandia.InterSpec/` |

So on macOS, for example, the full path to drop templates into is

```
~/Library/Application Support/sandia.InterSpec/act_shield_fit_templates/
```

### 3.3 Filename rules

A file in the writable directory is picked up only if its name contains one of these
three case-insensitive tokens:

- `tmplt.`     → e.g. `my_report.tmplt.html`
- `template.`  → e.g. `my_report.template.txt`
- `inja.`      → e.g. `my_report.inja.html`

The file's extension determines how InterSpec renders it: `.html` is treated as HTML
(rendered in a browser-style window), anything else is treated as plain text. The
"name" shown to the user in the dropdown is the filename with `tmplt` / `template`
stripped out (`inja` is *not* stripped) — so `my_report.tmplt.html` shows as
`my_report.html`, but `my_report.inja.html` shows as `my_report.inja.html`.

Recommended extensions: `.tmplt.html`, `.tmplt.txt`, `.tmplt.csv`, `.tmplt.md`.

Output file names are derived from the input spectrum file's leaf name.  When an input file
holds more than one candidate foreground spectrum and `--multi-sample-handling` is used to split
or sum it (see §3.4), an infix is added just before the extension so the results of each analysis
land in their own files: `foo.n42` becomes `foo_sample3.n42` for each sample, or `foo_summed.n42`
when the samples are summed.

### 3.4 Batch CLI overrides

The batch CLI (`InterSpec --batch-act-fit ...` or `--batch-peak-fit ...`) accepts:

- `--multi-sample-handling <auto|each|sum>` — how to treat an input file that holds more than
  one candidate foreground spectrum (i.e., more than one sample number not marked background,
  intrinsic activity, or calibration).  `auto` (the default) is the historical behaviour, and
  skips such a file with a warning.  `each` analyzes every candidate foreground sample
  separately, as if each was its own input file.  `sum` sums them into a single spectrum.
  Passthrough / search-mode files are never affected.
- `--file-report-template <name-or-path>` — a per-file template. May be given multiple
  times to render multiple per-file outputs. Accepted values: a template filename, an
  absolute path, or one of the built-in shorthands registered by
  `get_default_inja_env()` (e.g. `default-act-fit-html-results`,
  `default-peak-fit-html-results`).
- `--summary-report-template <name-or-path>` — a multi-file summary template (uses the
  `Files[]` wrapper described in §5.16 / §6.3). Same accepted values.
- `--report-template-include-dir <dir>` — an additional directory to search and to use
  as the search root for any `{% include %}` directives the template uses. The special
  tokens `default` (use the bundled directory) and `none` (disable file includes) are
  also accepted.

When run from the GUI ("Calculation Log → Custom Template…"), InterSpec auto-discovers
all valid templates in §3.2 and lists them in the dropdown.

### 3.5 Finding the writable directory from inside InterSpec

The easiest way to find — and open — the writable directory for your installation is
from InterSpec itself:

**Help → About InterSpec → Data** shows the resolved `<writableDataDirectory>` path
for your platform and build, with a button/link to open it in the system file browser
(Finder / Explorer / your DE's file manager). Create the `act_shield_fit_templates/`
sub-directory there if it doesn't already exist, and drop your custom template files
into it.

### 3.6 Editing workflow

To customise a default template:

1. Copy the bundled file (e.g. `act_fit.tmplt.html`) to the writable directory.

   On **macOS**:
   ```sh
   mkdir -p ~/Library/Application\ Support/sandia.InterSpec/act_shield_fit_templates
   cp /Applications/InterSpec.app/Contents/Resources/InterSpec_resources/static_text/ShieldSourceFitLog/act_fit.tmplt.html \
      ~/Library/Application\ Support/sandia.InterSpec/act_shield_fit_templates/my_act_fit.tmplt.html
   ```

   On **Windows** (PowerShell — adjust the install path if InterSpec lives elsewhere):
   ```powershell
   $dst = "$env:APPDATA\InterSpec\act_shield_fit_templates"
   New-Item -ItemType Directory -Force -Path $dst | Out-Null
   Copy-Item "C:\Program Files\InterSpec\resources\app\InterSpec_resources\static_text\ShieldSourceFitLog\act_fit.tmplt.html" `
             "$dst\my_act_fit.tmplt.html"
   ```

   If you're not sure where InterSpec is installed, use **Help → About InterSpec → Data**
   (see §3.5) to open the writable directory directly, then drop your edited copy in
   alongside the existing templates (you can also drag-and-drop a bundled template from
   InterSpec's installation directory into the writable directory using your file
   manager).

2. Edit `my_act_fit.tmplt.html` in your editor of choice.
3. Restart InterSpec (or reopen the Activity/Shielding tool) — the new template appears
   in the Calculation Log dropdown as `my_act_fit.html`.

## 4. How a report is rendered

The C++ side ([`src/BatchInfoLog.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/BatchInfoLog.cpp),
[`src/ShieldingSourceDisplay.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/ShieldingSourceDisplay.cpp),
[`src/BatchPeak.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/BatchPeak.cpp))
is a three-step pipeline that each caller wires up:

```cpp
// 1. Build the JSON payload (full data; see §5 or §6 for shape).
nlohmann::json data;
BatchInfoLog::shield_src_fit_results_to_json( results, drf, options, useBq, distance,
                                              geometry, data );
// (or, for peak-fit:)
BatchInfoLog::add_peak_fit_results_to_json( data, peak_fit_results );

// 2. Build an inja::Environment with the bundled templates and callbacks pre-registered.
//    `include_dir` is the search root for any {% include %} directives in user templates.
inja::Environment env = BatchInfoLog::get_default_inja_env( include_dir );

// 3. Render the chosen template against `data` using `env`.
env.write( "my_template.tmplt.html", data, out_stream );
```

`get_default_inja_env(include_dir)` enables `set_trim_blocks(true)`, registers the two
custom callbacks (§9), and pre-registers the eight bundled templates as named includes:

- `default-act-fit-html-results`     (= `act_fit.tmplt.html`)
- `default-act-fit-txt-results`      (= `std_fit_log.tmplt.txt`)
- `default-act-fit-html-summary`     (= `std_summary.tmplt.html`)
- `default-act-fit-csv-summary`      (= `std_summary.tmplt.csv`)
- `default-peak-fit-html-results`    (= `peak_fit.tmplt.html`)
- `default-peak-fit-txt-results`     (= `std_peak_fit_log.tmplt.txt`)
- `default-peak-fit-html-summary`    (= `std_peak_fit_summary.tmplt.html`)
- `default-peak-fit-csv-summary`     (= `std_peak_fit_summary.tmplt.csv`)

These can be invoked from a custom template via `{% include "default-..." %}`.

**Single-spectrum vs batch summary.** In single-spectrum mode the JSON described in §5 /
§6 is at the *top level* of `data`. In batch-summary mode the same payload is duplicated
inside each element of a top-level `Files[]` array, and `data` itself carries the
aggregate metadata (`InputFiles[]`, `AnalysisTime`, chart assets, etc.). Templates that
should work in both modes typically test for `exists("Files")` and dispatch.

## 5. JSON data model — Activity / Shielding fit

This section describes the JSON object built by
`BatchInfoLog::shield_src_fit_results_to_json()`. Every field below is referenced
directly by C++ code in `src/BatchInfoLog.cpp`; if you do not see a field in this section,
it does not exist in the payload.

### 5.1 Top-level keys

```json
{
  "Filepath": "/data/spectra/cs137_pb_shielded.n42",
  "Filename": "cs137_pb_shielded.n42",
  "ParentDir": "/data/spectra",
  "Success": true,
  "HasFitResults": true,
  "HasErrorMessage": false,
  "ErrorMessage": "",
  "HasWarnings": false,
  "Warnings": [],
  "FitChi2": 12.34,
  "EstimatedDistanceToMinimum": 0.0017,
  "NumberFcnCalls": 187,
  "NumDof": 7,
  "NumSources": 1,
  "HasTraceSource": false,
  "HasSelfAttenSource": false,
  "HasVolumetricSource": false,
  "AnySourceAgeFit": false,
  "AnyShieldingFit": true,
  "FixedGeometryDetector": false,
  "RawFitParameter": { "Values": [...], "Errors": [...] },
  "Detector": { ... },              // §5.2
  "ActShieldFitSetup": { ... },     // §5.3
  "foreground": { ... },            // §5.4
  "background": { ... },            // §5.4 (only if a background is loaded)
  "Sources": [ ... ],               // §5.6
  "Shieldings": { ... },            // §5.9
  "PullTrend": { ... },             // §5.18 (automatic diagnosis of the per-peak residual
                                    //  trend vs energy - "too much shielding", etc.)
  "VolumetricEff": { ... },         // §5.19 (which per-element efficiency the volumetric-
                                    //  source integration actually used, and why)
  "PeaksUsedForActivityFitting": { "Peaks": [ ... ] },  // §5.12
  "PeakToModelComparison": { "UsedPeaks": [ ... ] },    // §5.13
  "AnyNotFitPeakMda": true,
  "NotFitPeaks": { ... },                               // §6.4 (exemplar peaks that did not
                                                        //  converge, with their detection
                                                        //  limits - §6.6.1)
  "SupplementalPeakInfo": { ... },                      // §5.17 (per-peak detection-limit checks,
                                                        //  and the activities implied by peaks
                                                        //  not used in the fit)
  "D3_JS": "<minified d3.v3.js>",                       // §5.14
  "SpectrumChart_JS": "...",
  "SpectrumChart_CSS": "...",
  "ShieldingSourceFitPlot_JS": "...",
  "ShieldingSourceFitPlot_CSS": "...",
  "InterSpecCompileDate": "May 20 2026",                // §5.15
  "InterSpecCompileDateIso": "2026-05-20",
  "AnalysisTime": "2026-05-20T14:30:45",
  "CurrentWorkingDirectory": "/home/wcjohns/work",
  "InterSpecExecutablePath": "/Applications/InterSpec.app/Contents/MacOS/InterSpec"
}
```

| Field | Type | Notes |
|---|---|---|
| `Filepath`                     | string  | Full path to the input spectrum file. |
| `Filename`                     | string  | Leaf name identifying this analysis, and matching the output files written for it. Equals the input file's leaf name, unless multi-sample handling split or summed the file, in which case it carries the same `_sampleN` / `_summed` infix the output files do. |
| `SourceFilename`               | string  | The input file's leaf name, never decorated. |
| `AnalysisLabel`                | string  | Human-readable label for this analysis, e.g. `foo.n42 (sample 3)`. |
| `IsSplitFromMultiSampleFile`   | bool    | True when this analysis is one of several taken from a single input file. |
| `ForegroundSampleNumbers`      | array[int] | Sample numbers used as the foreground. Only present when they were explicitly determined; more than one means they were summed. |
| `ParentDir`                    | string  | Directory containing the input file. |
| `Success`                      | bool    | Fit reached `Final` status (i.e. converged and not cancelled / aborted). |
| `HasFitResults`                | bool    | At least some fit results exist (even if `Success=false`). |
| `HasErrorMessage`              | bool    | `ErrorMessage` is non-empty. Set by the GUI / batch driver wrapper. |
| `ErrorMessage`                 | string  | Concatenated error text from the fit. Only present when there *is* an error message — guard with `HasErrorMessage` or `exists("ErrorMessage")`. |
| `ErrorMessages`                | array[string] | Per-error array emitted by `shield_src_fit_results_to_json` itself. Only present if non-empty; the singular `ErrorMessage` above (from the wrapping caller) is the more reliable field for templates. |
| `HasWarnings`                  | bool    | `Warnings` is non-empty. Always emitted. |
| `Warnings`                     | array[string] | Human-readable warning strings. Always emitted (may be empty). |
| `FitChi2`                      | number  | Final χ² value reported by Minuit2. |
| `EstimatedDistanceToMinimum`   | number  | EDM from Minuit2; small means well-converged. |
| `NumberFcnCalls`               | int     | Number of objective-function evaluations. |
| `NumDof`                       | int     | Degrees of freedom (number of free parameters). |
| `NumSources`                   | int     | Number of distinct nuclides being fit. |
| `HasTraceSource`               | bool    | At least one source is a trace source (§5.6). |
| `HasSelfAttenSource`           | bool    | At least one source is self-attenuating (§5.6). |
| `HasVolumetricSource`          | bool    | `HasTraceSource OR HasSelfAttenSource`. |
| `AnySourceAgeFit`              | bool    | At least one source's age was fitted. |
| `AnyShieldingFit`              | bool    | At least one shielding dimension or property was fitted. |
| `FixedGeometryDetector`        | bool    | The DRF has a fixed-geometry efficiency (no `Distance`/`Geometry`). |
| `RawFitParameter.Values`       | array[number] | Raw parameter vector from Minuit2 (implementation-defined ordering — for diagnostics). |
| `RawFitParameter.Errors`       | array[number] | Raw parameter uncertainties (paired with `.Values`). |

### 5.2 `Detector`

```json
"Detector": {
  "Name": "HPGe-Detective",
  "Description": "Detective-X HPGe",
  "Diameter": "7.62 cm", "Diameter_mm": 76.2, "Diameter_cm": 7.62, "Diameter_m": 0.0762, "Diameter_inch": 3.0,
  "Radius":   "3.81 cm", "Radius_mm":   38.1, "Radius_cm":   3.81, "Radius_m":   0.0381, "Radius_inch":   1.5,
  "Setback":   "0 mm",   "Setback_cm":  0.0,
  "FixedGeometry": false
}
```

| Field | Type | Units | Notes |
|---|---|---|---|
| `Name`              | string  |     | Detector name as stored in the DRF. |
| `Description`       | string  |     | Free-text description from the DRF. |
| `Diameter`          | string  |     | Human-readable, "best units" formatting. |
| `Diameter_mm`       | number  | mm  | |
| `Diameter_cm`       | number  | cm  | |
| `Diameter_m`        | number  | m   | |
| `Diameter_inch`     | number  | in  | |
| `Radius`            | string  |     | Human-readable, "best units" formatting. |
| `Radius_mm`         | number  | mm  | |
| `Radius_cm`         | number  | cm  | |
| `Radius_m`          | number  | m   | |
| `Radius_inch`       | number  | in  | |
| `Setback`           | string  |     | Distance behind the detector face. Only present when setback > 0. |
| `Setback_cm`        | number  | cm  | Only present when setback > 0. |
| `FixedGeometry`     | bool    |     | True if this DRF is a fixed-geometry efficiency. |

### 5.3 `ActShieldFitSetup`

```json
"ActShieldFitSetup": {
  "Distance":      "1.00 m",
  "Distance_mm":   1000.0,
  "Distance_cm":   100.0,
  "Distance_m":    1.0,
  "Distance_km":   0.001,
  "Distance_inch": 39.37,
  "Distance_feet": 3.281,
  "Geometry": "Spherical",
  "FixedGeometryType": "",
  "FitOptions": {
    "InterferenceCorrection": true,
    "AttenuateForAir":        true,
    "DecayDuringMeasurement": false,
    "MultithreadSelfAttenCalc": true,
    "PhotopeakClusterSigma":  1.5,
    "BackgroundPeakSubtract": false,
    "ElementNuclidesSameAge": true
  }
}
```

| Field | Type | Units | Notes |
|---|---|---|---|
| `Distance`        | string  |     | Best-units distance. Only present when not fixed-geometry. |
| `Distance_mm`     | number  | mm  | Only present when not fixed-geometry. |
| `Distance_cm`     | number  | cm  | Only present when not fixed-geometry. |
| `Distance_m`      | number  | m   | Only present when not fixed-geometry. |
| `Distance_km`     | number  | km  | Only present when not fixed-geometry. |
| `Distance_inch`   | number  | in  | Only present when not fixed-geometry. |
| `Distance_feet`   | number  | ft  | Only present when not fixed-geometry. |
| `Geometry`        | string  |     | `"Spherical"`, `"CylinderEndOn"`, `"CylinderSideOn"`, or `"Rectangular"`. Only present when not fixed-geometry. |
| `FixedGeometryType` | string |    | Human-readable description (`"total activity"`, `"activity per square centimeter"`, `"activity per square meter"`, `"activity per gram"`). Only present when fixed-geometry. |
| `FitOptions.InterferenceCorrection`   | bool   |  | Multiple nuclides may contribute to the same photopeak. |
| `FitOptions.AttenuateForAir`          | bool   |  | Air attenuation included between source and detector. |
| `FitOptions.DecayDuringMeasurement`   | bool   |  | Branching-ratios corrected for decay during the measurement. |
| `FitOptions.MultithreadSelfAttenCalc` | bool   |  | Self-attenuation integrals are parallelised. |
| `FitOptions.PhotopeakClusterSigma`    | number |  | Energies within this many σ are clustered into a single photopeak. |
| `FitOptions.BackgroundPeakSubtract`   | bool   |  | Background peaks are subtracted from foreground peaks before fitting. |
| `FitOptions.ElementNuclidesSameAge`   | bool   |  | Different isotopes of the same element share one age parameter. |

### 5.4 `foreground` and `background`

Both objects have an identical structure (built by `add_hist_to_json`). `background` is
only present if a background spectrum was loaded; `foreground` is always present when
`HasFitResults=true`.

```json
"foreground": {
  "LiveTime": "300.00 s",       "LiveTime_s": 300.0,
  "RealTime": "301.20 s",       "RealTime_s": 301.2,
  "DeadTime": "1.20 s",         "DeadTime_s": 1.2,    "DeadTime_percent": 0.398,
  "StartTime":     "2026-05-20T14:30:45.123456",
  "StartTime_iso": "2026-05-20T14:30:45",
  "StartTime_vax": "20-MAY-2026 14:30:45.12",
  "StartTimeIsValid": true,
  "LowerSpectrumEnergy": 0.0,   "UpperSpectrumEnergy": 3000.0,
  "NumberChannels": 8192,
  "GammaSum": 1.23e6,           "GammaCps": 4100.0,
  "SampleNumbers": [1],
  "Filename": "cs137_pb_shielded.n42",
  "spectrum": { ... },           // raw D3 spectrum JSON (only if SpecUtils built with D3 support)
  "HasGps": true,
  "Longitude": -106.5,
  "Latitude":   35.7,
  "HasNeutrons": false,
  "NeutronCounts": 0.0, "NeutronLiveTime": 0.0, "NeutronCps": 0.0,
  "SpectrumTitle":        "",
  "DoseRate_uSvPerHour":  null,
  "ExposureRate_mRPerHour": null,
  "DetectorTypeDesc": "HPGe",
  "DetectorName":     "Detective-X 1234",
  "SpectrumRemarks":  "",
  "SpectrumParseWarnings": "",
  "IsDerivedData":         false,
  "DerivedIoiSum":         false,
  "DerivedUsedForRidAnalysis": false,
  "DerivedProcessedFurther":   false,
  "DerivedBackgroundSub":      false,
  "DerivedIsBackground":       false,
  "EnergyCal": { ... },          // §5.5
  "InstrumentModel": "...",
  "SerialNumber":    "...",
  "Manufacturer":    "...",
  "InstrumentType":  "...",
  "DetectorType":    "...",
  "NumberRecordsInFile": 1,
  "RemarksInFile":       "...",
  "ParseWarningsForFile":"...",
  "HasInstrumentRid":    false,
  "InstrumentRidSummary":"",
  "Normalization": 1.0           // background only; foreground does not have this field
}
```

| Field | Type | Units | Notes |
|---|---|---|---|
| `LiveTime`            | string  |     | Best-units formatting (e.g. `"123.45 s"`). |
| `LiveTime_s`          | number  | s   | |
| `RealTime`            | string  |     | |
| `RealTime_s`          | number  | s   | |
| `DeadTime`            | string  |     | `RealTime − LiveTime`, best-units. |
| `DeadTime_s`          | number  | s   | |
| `DeadTime_percent`    | number  | %   | `100·DeadTime / RealTime`. |
| `StartTime`           | string  |     | ISO-8601 extended (with microseconds). |
| `StartTime_iso`       | string  |     | ISO-8601 basic (seconds resolution). |
| `StartTime_vax`       | string  |     | VAX-style timestamp. |
| `StartTimeIsValid`    | bool    |     | False if the start-time field was missing from the input file. |
| `LowerSpectrumEnergy` | number  | keV | Lower edge of the lowest channel. |
| `UpperSpectrumEnergy` | number  | keV | Upper edge of the highest channel. |
| `NumberChannels`      | int     |     | |
| `GammaSum`            | number  | cnt | Total gamma counts in the spectrum. |
| `GammaCps`            | number  | s⁻¹ | `GammaSum / LiveTime_s`. |
| `SampleNumbers`       | array[int] |  | Sample numbers used (often `[1]`). |
| `Filename`            | string  |     | Source filename for this spectrum. |
| `spectrum`            | object  |     | Raw D3 chart JSON, suitable to embed in a `<script>` block. Only present when InterSpec is built with D3 chart support. |
| `HasGps`              | bool    |     | |
| `Longitude`           | number  | °   | Only present if `HasGps`. |
| `Latitude`            | number  | °   | Only present if `HasGps`. |
| `HasNeutrons`         | bool    |     | |
| `NeutronCounts`       | number  | cnt | Only present if `HasNeutrons`. |
| `NeutronLiveTime`     | number  | s   | Only present if `HasNeutrons`. |
| `NeutronCps`          | number  | s⁻¹ | Only present if `HasNeutrons`. |
| `SpectrumTitle`       | string  |     | Only present if non-empty in the file. |
| `DoseRate_uSvPerHour` | number  | µSv/h | Only present if reported by the instrument. |
| `ExposureRate_mRPerHour` | number | mR/h | Only present if reported by the instrument. |
| `DetectorTypeDesc`    | string  |     | Detector type description from the file (e.g. `"HPGe"`). |
| `DetectorName`        | string  |     | Detector model / serial / name from the file. |
| `SpectrumRemarks`     | string  |     | Per-record remarks. |
| `SpectrumParseWarnings` | string |    | Per-record parse warnings. |
| `IsDerivedData`        | bool   |     | Spectrum is a derived (processed) record, not raw measurement. |
| `DerivedIoiSum`        | bool   |     | Derived item-of-interest sum. |
| `DerivedUsedForRidAnalysis` | bool |  | Derived record was used for the instrument's built-in RID. |
| `DerivedProcessedFurther`   | bool |  | Derived record was further processed. |
| `DerivedBackgroundSub`      | bool |  | Derived record is background-subtracted. |
| `DerivedIsBackground`       | bool |  | Derived record is itself the background. |
| `EnergyCal`           | object  |     | See §5.5. |
| `InstrumentModel`     | string  |     | Always emitted when the spectrum's source file metadata is available (may be an empty string). |
| `SerialNumber`        | string  |     | Always emitted when available; may be empty. |
| `Manufacturer`        | string  |     | Always emitted when available; may be empty. |
| `InstrumentType`      | string  |     | Always emitted when available; may be empty. |
| `DetectorType`        | string  |     | Always emitted when available; may be empty. |
| `NumberRecordsInFile` | int     |     | Always emitted when file metadata is available. |
| `RemarksInFile`       | string  |     | File-level remarks. Always emitted when file metadata is available; may be empty. |
| `ParseWarningsForFile`| string  |     | File-level parse warnings. Always emitted when file metadata is available; may be empty. |
| `HasInstrumentRid`    | bool    |     | The file embedded a manufacturer-RID result. Only present if the file had RID metadata. |
| `InstrumentRidSummary`| string  |     | Human-readable summary of the manufacturer RID. Only present if `HasInstrumentRid`. |
| `Normalization`       | number  |     | **Background only.** Live-time normalization factor applied when subtracting background from foreground. |

### 5.5 `foreground.EnergyCal` (and `background.EnergyCal`)

```json
"EnergyCal": {
  "NumChannels": 8192,
  "Type": "Polynomial",
  "LowerEnergy": 0.0,
  "UpperEnergy": 3000.0,
  "Coefficients": [0.0, 0.366, 1.2e-7],
  "DeviationPairs": [[122.06, 0.0], [661.66, 0.31], [1332.5, 0.18]]
}
```

| Field | Type | Notes |
|---|---|---|
| `NumChannels`     | int           | Channel count. |
| `Type`            | string        | `"Polynomial"`, `"FullRangeFraction"`, `"LowerChannelEdge"`, or `"Invalid"`. |
| `LowerEnergy`     | number (keV)  | Lowest calibrated energy. Only present if `Type != "Invalid"`. |
| `UpperEnergy`     | number (keV)  | Highest calibrated energy. Only present if `Type != "Invalid"`. |
| `Coefficients`    | array[number] | Polynomial or FRF coefficients. Only present for `"Polynomial"` / `"FullRangeFraction"`. |
| `DeviationPairs`  | array[[number,number]] | Pairs of `[energy_keV, deviation_keV]`. Only present if any are defined. |

For `"LowerChannelEdge"` calibration, the per-channel lower edges live inside `Coefficients`.

### 5.6 `Sources[]`

Each element of the top-level `Sources[]` array describes one fitted nuclide. The
structure is built by `BatchInfoLog::add_basic_src_details()` plus a per-source
`PeaksThisNucContributesTo[]` array.

```json
{
  "Nuclide": "Cs-137",
  "Activity":     "12.34 µCi",
  "Activity_bq":  4.566e5,  "Activity_kBq": 456.6, "Activity_MBq": 0.4566, "Activity_GBq": 4.566e-4,
  "Activity_ci":  1.234e-5, "Activity_mCi": 0.01234, "Activity_uCi": 12.34, "Activity_pCi": 1.234e7,
  "ActivityPostFix": "",
  "ActivityIsFit":   true,
  "ActivityUncert":  "0.5 µCi",
  "ActivityUncert_bq":  1.85e4, "ActivityUncert_kBq": 18.5, "ActivityUncert_MBq": 0.0185, "ActivityUncert_GBq": 1.85e-5,
  "ActivityUncert_ci":  5.0e-7, "ActivityUncert_mCi": 5.0e-4, "ActivityUncert_uCi": 0.5, "ActivityUncert_pCi": 5.0e5,
  "ActivityUncertPercent": "4.05",
  "NuclideMass": "143 ng",
  "Age": "30.0 y",
  "AgeSeconds": 9.467e8, "AgeDays": 1.096e4, "AgeYears": 30.0,
  "AgeIsFittable": false,
  "AgeIsFit": false,
  "AgeDefiningNuclide": "Cs-137",
  "IsTraceSource": false,
  "IsSelfAttenSource": false,
  "PeaksThisNucContributesTo": [ ... ]      // §5.7
}
```

| Field | Type | Units | Notes |
|---|---|---|---|
| `Nuclide`                  | string  |        | E.g. `"Cs-137"`, `"Co-60"`. |
| `Activity`                 | string  |        | Best-units string, with `ActivityPostFix` appended. |
| `Activity_bq`              | number  | Bq     | |
| `Activity_kBq`             | number  | kBq    | |
| `Activity_MBq`             | number  | MBq    | |
| `Activity_GBq`             | number  | GBq    | |
| `Activity_ci`              | number  | Ci     | |
| `Activity_mCi`             | number  | mCi    | |
| `Activity_uCi`             | number  | µCi    | Yes, the *field* name is `uCi` but the *unit* is µCi (microcuries). |
| `Activity_pCi`             | number  | pCi    | |
| `ActivityPostFix`          | string  |        | `""` for point sources; `"/cm^2"` / `"/g"` / `"/cm^3"` / `"/m2"` for fixed-geometry detectors. |
| `ActivityIsFit`            | bool    |        | True if activity was a free parameter. |
| `ActivityUncert`           | string  |        | Best-units uncertainty string. Only present if `ActivityIsFit`. |
| `ActivityUncert_bq`        | number  | Bq     | Only present if `ActivityIsFit`. |
| `ActivityUncert_kBq`       | number  | kBq    | Only present if `ActivityIsFit`. |
| `ActivityUncert_MBq`       | number  | MBq    | Only present if `ActivityIsFit`. |
| `ActivityUncert_GBq`       | number  | GBq    | Only present if `ActivityIsFit`. |
| `ActivityUncert_ci`        | number  | Ci     | Only present if `ActivityIsFit`. |
| `ActivityUncert_mCi`       | number  | mCi    | Only present if `ActivityIsFit`. |
| `ActivityUncert_uCi`       | number  | µCi    | Only present if `ActivityIsFit`. |
| `ActivityUncert_pCi`       | number  | pCi    | Only present if `ActivityIsFit`. |
| `ActivityUncertPercent`    | string  | %      | `100 · ActivityUncert / Activity`, formatted via `printCompact(...,4)`. Only present if `ActivityIsFit`. |
| `NuclideMass`              | string  |        | Best-units mass (ng / µg / mg / g …). |
| `Age`                      | string  |        | Best-units age. |
| `AgeSeconds`               | number  | s      | |
| `AgeDays`                  | number  | d      | |
| `AgeYears`                 | number  | y      | |
| `AgeIsFittable`            | bool    |        | True if the nuclide's age is meaningfully variable (has decay daughters whose branching-ratios change with age). |
| `AgeIsFit`                 | bool    |        | True if age was a free parameter. |
| `AgeUncert`                | string  |        | Only present if `AgeIsFit`. |
| `AgeUncertSeconds`         | number  | s      | Only present if `AgeIsFit`. |
| `AgeUncertDays`            | number  | d      | Only present if `AgeIsFit`. |
| `AgeUncertYears`           | number  | y      | Only present if `AgeIsFit`. |
| `AgeDefiningNuclide`       | string  |        | Symbol of the age-defining nuclide for this element (when `ElementNuclidesSameAge` couples ages). Only present if non-null. |
| `IsTraceSource`            | bool    |        | True if this is a trace source (see below). |
| `TraceActivityType`        | string  |        | `"TotalActivity"`, `"ActivityPerCm3"`, `"ActivityPerGram"`, or `"ExponentialDistribution"`. Only present if `IsTraceSource`. |
| `TraceDisplayActivity`     | string  |        | Best-units trace-display activity with postfix. Only present if `IsTraceSource`. |
| `TraceDisplayActivity_bq`  | number  | Bq     | Only present if `IsTraceSource`. |
| `TraceDisplayActivity_kBq` | number  | kBq    | Only present if `IsTraceSource`. |
| `TraceDisplayActivity_MBq` | number  | MBq    | Only present if `IsTraceSource`. |
| `TraceDisplayActivity_GBq` | number  | GBq    | Only present if `IsTraceSource`. |
| `TraceDisplayActivity_ci`  | number  | Ci     | Only present if `IsTraceSource`. |
| `TraceDisplayActivity_mCi` | number  | mCi    | Only present if `IsTraceSource`. |
| `TraceDisplayActivity_uCi` | number  | µCi    | Only present if `IsTraceSource`. |
| `TraceDisplayActivity_pCi` | number  | pCi    | Only present if `IsTraceSource`. |
| `TraceActivityPostFix`     | string  |        | `""`, `"/cm^3"`, `"/g"`, or `"/m2 exp"`. Only present if `IsTraceSource`. |
| `TraceDisplayActivityUncert`     | string |    | Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceDisplayActivityUncert_bq`  | number | Bq | Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceDisplayActivityUncert_kBq` | number | kBq| Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceDisplayActivityUncert_MBq` | number | MBq| Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceDisplayActivityUncert_GBq` | number | GBq| Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceDisplayActivityUncert_ci`  | number | Ci | Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceDisplayActivityUncert_mCi` | number | mCi| Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceDisplayActivityUncert_uCi` | number | µCi| Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceDisplayActivityUncert_pCi` | number | pCi| Only present if `IsTraceSource` and `ActivityIsFit`. |
| `TraceRelaxationLength`         | string |     | Only present if `TraceActivityType == "ExponentialDistribution"`. |
| `TraceRelaxationLength_mm`      | number | mm  | Only present if `TraceActivityType == "ExponentialDistribution"`. |
| `TraceRelaxationLength_cm`      | number | cm  | Only present if `TraceActivityType == "ExponentialDistribution"`. |
| `TraceRelaxationLength_m`       | number | m   | Only present if `TraceActivityType == "ExponentialDistribution"`. |
| `TraceRelaxationLength_inch`    | number | in  | Only present if `TraceActivityType == "ExponentialDistribution"`. |
| `IsSelfAttenSource`             | bool   |     | True if this nuclide lives inside one of the shielding layers. |
| `SelfAttenShieldIndex`          | int    |     | Index into `Shieldings.Shields[]` of the host layer. Only present if `IsSelfAttenSource`. |
| `SelfAttenShieldName`           | string |     | Name of the host layer. Only present if `IsSelfAttenSource`. |
| `SelfAttenIsVariableMassFrac`   | bool   |     | The mass fraction of this nuclide in the host layer was a free parameter. Only present if `IsSelfAttenSource`. |
| `SelfAttenMassFrac`             | number | 0–1 | Mass fraction in the host layer. Only present if `IsSelfAttenSource`. |
| `SelfAttenMassFracUncert`       | number | 0–1 | Mass-fraction uncertainty. Only present if `IsSelfAttenSource` and `SelfAttenIsVariableMassFrac`. |
| `PeaksThisNucContributesTo`     | array[object] | | See §5.7. May be empty if no peaks were assigned to this nuclide. |

### 5.7 `Sources[].PeaksThisNucContributesTo[]`

Each element describes one fitted peak that the source contributes gammas to. Built by
`BatchInfoLog::add_basic_peak_info()` plus a nested `ThisNucsGammasForPeak[]` array
(§5.8).

| Field | Type | Units | Notes |
|---|---|---|---|
| `Energy`                     | string | keV   | Formatted with `%.2f` (2 decimal places, e.g. `"661.66"`). |
| `Energy_keV`                 | number | keV   | Peak centroid. |
| `DecayParticleEnergy`        | string | keV   | Formatted with `%.3f` (3 decimal places). |
| `DecayParticleEnergy_keV`    | number | keV   | |
| `AssignedNuclide`            | string |       | Symbol of the nuclide assigned to the peak (may differ from the outer source if the peak is shared). |
| `FWHM`                       | number | keV   | Full width at half-maximum. |
| `Counts`                     | number | cnt   | Fitted peak area. |
| `CountsStr`                  | string |       | `printCompact(Counts, 4)`. |
| `CountsUncert`               | number | cnt   | |
| `CountsUncertStr`            | string |       | |
| `Cps`                        | number | s⁻¹   | `Counts / foreground.LiveTime_s`. |
| `CpsStr`                     | string |       | |
| `CpsUncert`                  | number | s⁻¹   | |
| `BackgroundCounts`           | number | cnt   | Background-spectrum contribution at this peak. Only present if non-zero. |
| `BackgroundCountsStr`        | string |       | Only present if non-zero. |
| `BackgroundCountsUncert`     | number | cnt   | Only present if non-zero and uncertainty available. |
| `BackgroundCountsUncertStr`  | string |       | Only present if non-zero and uncertainty available. |
| `SignalCounts`               | number | cnt   | `Counts − BackgroundCounts`. |
| `SignalCountsStr`            | string |       | |
| `SignalCountsUncert`         | number | cnt   | |
| `SignalCountsUncertStr`      | string |       | |
| `PredictedCounts`            | number | cnt   | What the fitted model predicts for this peak. |
| `PredictedNumSigmaOff`       | number | σ     | `(SignalCounts − PredictedCounts) / SignalCountsUncert`. |
| `ObservedOverPredicted`      | number |       | `SignalCounts / PredictedCounts`. |
| `ObservedOverPredictedUncert`| number |       | Propagated uncertainty on the ratio. |
| `DetectorSolidAngleFraction` | number | 0–1   | Fraction of 4π subtended by the detector. |
| `DetectorIntrinsicEff`       | number | 0–1   | Intrinsic photopeak efficiency. |
| `DetectorEff`                | number | 0–1   | `DetectorSolidAngleFraction × DetectorIntrinsicEff`. |
| `ShieldAttenuations`         | array[number] | | One factor per shielding layer (`Shieldings.Shields[]` order). |
| `AttenuationByShieldingFactor` | number |     | Product of `ShieldAttenuations[]`. |
| `AttenuationByAirFactor`     | number |       | Air attenuation factor between source and detector. |
| `AttenuationTotalFactor`     | number |       | `AttenuationByShieldingFactor × AttenuationByAirFactor`. |
| `ThisNucsGammasForPeak`      | array[object] | | See §5.8. |

### 5.8 `Sources[].PeaksThisNucContributesTo[].ThisNucsGammasForPeak[]`

Each element describes one discrete gamma line from the outer source that contributes to
the peak. Built by `BatchInfoLog::add_gamma_info_for_peak()`. Note: this object **also
contains all of the `add_basic_src_details()` fields** for the outer source (i.e. every
`Activity_*`, `Age*`, trace/self-atten field listed in §5.6 is duplicated here for
convenience), so you can write a single loop over gammas without joining back to
`Sources[]`.

Gamma-specific fields:

| Field | Type | Units | Notes |
|---|---|---|---|
| `Nuclide`              | string |       | Symbol of the gamma's parent nuclide. |
| `Energy`               | string | keV   | Formatted with `%.3f` (3 decimal places). |
| `Energy_keV`           | number | keV   | |
| `BranchingRatio`       | number | 0–1   | Branching ratio at the source's age. |
| `BranchingRatioStr`    | string |       | |
| `PredictedCounts`      | number | cnt   | Model prediction *for this gamma's contribution* to the parent peak. |
| `PredictedCountsStr`   | string |       | |
| `SourcePhotonsCps`     | number | γ/s   | Photons per second emitted at the source for this gamma. |
| `SourcePhotonsCpsStr`  | string |       | |
| `SourcePhotons`        | number | γ     | `SourcePhotonsCps · foreground.LiveTime_s`. |
| `SourcePhotonsStr`     | string |       | |
| `HasDecayCorrection`   | bool   |       | True if a decay-during-measurement correction was applied. |
| `DecayCorrection`      | number |       | Correction factor (`(1 − exp(−λT)) / (λT)`-like). Only present if `HasDecayCorrection`. |
| `DecayCorrectionStr`   | string |       | Only present if `HasDecayCorrection`. |

### 5.9 `Shieldings` and `Shieldings.Shields[]`

```json
"Shieldings": {
  "Geometry": "Spherical",
  "NumberShieldings": 2,
  "NumDimensions": 1,
  "DimensionMeanings": ["Radius"],
  "Shields": [ ... ]
}
```

| Field | Type | Notes |
|---|---|---|
| `Geometry`         | string         | `"Spherical"`, `"CylinderEndOn"`, `"CylinderSideOn"`, `"Rectangular"`. |
| `NumberShieldings` | int            | Length of `Shields[]`. |
| `NumDimensions`    | int            | 1 (sphere), 2 (cylinder), or 3 (rectangular). |
| `DimensionMeanings`| array[string]  | Per-axis labels, e.g. `["Radius"]`, `["Radius","Length"]`, `["Width","Height","Depth"]`. |
| `Shields`          | array[object]  | One element per layer, ordered innermost → outermost. |

Each element of `Shields[]`:

```json
{
  "Name": "Lead",
  "ShieldingNumber": 0,
  "IsGeneric": false,
  "Formula": "Pb",
  "Density_gPerCm3": 11.34,
  "DensityStr": "11.34 g/cm3",
  "Thickness":      "5.00 mm",
  "Thickness_mm":   5.0, "Thickness_cm": 0.5, "Thickness_m": 0.005, "Thickness_inch": 0.1969, "Thickness_feet": 0.01641,
  "VolumeCm3": 12.57, "VolumeUncertCm3": 0.21,
  "InnerRadius": "0 mm",  "OuterRadius": "5.0 mm",
  "InnerDims":   ["0 mm"], "OuterDims": ["5.0 mm"], "Thicknesses": ["5.0 mm"], "ThicknessesUncerts": ["0.08 mm"],
  "InnerDims_mm":         [0.0],   "OuterDims_mm":         [5.0],   "Thicknesses_mm":         [5.0],   "ThicknessesUncerts_mm":   [0.08],
  "InnerDims_cm":         [0.0],   "OuterDims_cm":         [0.5],   "Thicknesses_cm":         [0.5],   "ThicknessesUncerts_cm":   [0.008],
  "InnerDims_m":          [0.0],   "OuterDims_m":          [0.005], "Thicknesses_m":          [0.005], "ThicknessesUncerts_m":    [8e-5],
  "InnerDims_inch":       [0.0],   "OuterDims_inch":       [0.197], "Thicknesses_inch":       [0.197], "ThicknessesUncerts_inch": [0.003],
  "NumDimensions": 1,
  "DimensionIsFit": [true],
  "FitAnyMassFraction": false,
  "SelfAttenSources": [ ... ],  // §5.10
  "TraceSources":     [ ... ]   // §5.11
}
```

| Field | Type | Units | Notes |
|---|---|---|---|
| `Name`                  | string  |          | Material name (`"Lead"`, `"Steel"`, or a user-defined material). For generic shieldings, a generated label. |
| `ShieldingNumber`       | int     |          | 0-based index in the array. |
| `Geometry`              | string  |          | Per-shield geometry, same vocabulary as `Shieldings.Geometry`: `"Spherical"`, `"CylinderEndOn"`, `"CylinderSideOn"`, `"Rectangular"`. (Currently always matches the top-level `Shieldings.Geometry`.) |
| `IsGeneric`             | bool    |          | True for AN/AD generic shielding; false for material-defined shielding. |
| `Formula`               | string  |          | Chemical formula of the material. Only present if `!IsGeneric`. |
| `Density_gPerCm3`       | number  | g/cm³    | Only present if `!IsGeneric`. |
| `DensityStr`            | string  |          | Best-units density. Only present if `!IsGeneric`. |
| `AN`                    | string  |          | Effective atomic number, formatted. Only present if `IsGeneric` and AN > 0. |
| `AD`                    | string  |          | Areal density, formatted. Only present if `IsGeneric` and AD > 0. |
| `FitAN`                 | bool    |          | True if AN was fit. Only present if `IsGeneric`. |
| `FitAD`                 | bool    |          | True if AD was fit. Only present if `IsGeneric`. |
| `Thickness`             | string  |          | Best-units scalar thickness (geometry-appropriate). |
| `Thickness_mm`          | number  | mm       | |
| `Thickness_cm`          | number  | cm       | |
| `Thickness_m`           | number  | m        | |
| `Thickness_inch`        | number  | in       | |
| `Thickness_feet`        | number  | ft       | |
| `VolumeCm3`             | number  | cm³      | Volume of this layer. |
| `VolumeUncertCm3`       | number  | cm³      | |
| `InnerRadius`           | string  |          | Best-units inner radius (sphere/cylinder). |
| `OuterRadius`           | string  |          | Best-units outer radius (sphere/cylinder). |
| `InnerDims`             | array[string] |    | Per-axis inner edges, best-units strings (length = `NumDimensions`). |
| `OuterDims`             | array[string] |    | Per-axis outer edges, best-units strings. |
| `Thicknesses`           | array[string] |    | Per-axis thickness strings. |
| `ThicknessesUncerts`    | array[string] |    | Per-axis thickness uncertainty strings. |
| `InnerDims_mm`          | array[number] | mm | Per-axis inner edge in mm. |
| `InnerDims_cm`          | array[number] | cm | |
| `InnerDims_m`           | array[number] | m  | |
| `InnerDims_inch`        | array[number] | in | |
| `OuterDims_mm`          | array[number] | mm | |
| `OuterDims_cm`          | array[number] | cm | |
| `OuterDims_m`           | array[number] | m  | |
| `OuterDims_inch`        | array[number] | in | |
| `Thicknesses_mm`        | array[number] | mm | |
| `Thicknesses_cm`        | array[number] | cm | |
| `Thicknesses_m`         | array[number] | m  | |
| `Thicknesses_inch`      | array[number] | in | |
| `ThicknessesUncerts_mm` | array[number] | mm | |
| `ThicknessesUncerts_cm` | array[number] | cm | |
| `ThicknessesUncerts_m`  | array[number] | m  | |
| `ThicknessesUncerts_inch` | array[number] | in | |
| `NumDimensions`         | int     |          | 1, 2, or 3. |
| `DimensionIsFit`        | array[bool] |      | Per-axis: was the thickness along this axis a free parameter? |
| `FitAnyMassFraction`    | bool    |          | True if any embedded self-atten source had a fit mass fraction. |
| `SelfAttenSources`      | array[object] |    | Nuclides embedded in this layer. Only present if non-empty. See §5.10. |
| `TraceSources`          | array[object] |    | Trace sources hosted in this layer. Only present if non-empty. See §5.11. |

### 5.10 `Shieldings.Shields[].SelfAttenSources[]`

Each element duplicates **all `add_basic_src_details()` fields** from §5.6 (so every
`Activity_*`, `Age*`, etc. is available), plus the following self-atten specifics:

| Field | Type | Notes |
|---|---|---|
| `Nuclide`                 | string  | Symbol. |
| `IsFittingMassFraction`   | bool    | The mass fraction was a free parameter. |
| `MassFraction`            | number  | Mass fraction (0–1) of this nuclide in the host layer. |
| `MassFractionStr`         | string  | Formatted via `printFixed(..., 6)`. |
| `MassFractionUncert`      | number  | Only present if `IsFittingMassFraction`. |
| `MassFractionUncertStr`   | string  | Only present if `IsFittingMassFraction`. |

### 5.11 `Shieldings.Shields[].TraceSources[]`

Each element duplicates **all `add_basic_src_details()` fields** from §5.6, plus:

| Field | Type | Units | Notes |
|---|---|---|---|
| `Nuclide`                | string  |       | Symbol. |
| `TraceSourceType`        | string  |       | `"TotalActivity"`, `"ActivityPerCm3"`, `"ActivityPerGram"`, `"ExponentialDistribution"`. |
| `RelaxationLength`       | string  |       | Only present if `TraceSourceType == "ExponentialDistribution"`. |
| `RelaxationLength_mm`    | number  | mm    | Only present if `TraceSourceType == "ExponentialDistribution"`. |
| `RelaxationLength_cm`    | number  | cm    | Only present if `TraceSourceType == "ExponentialDistribution"`. |
| `RelaxationLength_m`     | number  | m     | Only present if `TraceSourceType == "ExponentialDistribution"`. |
| `RelaxationLength_inch`  | number  | in    | Only present if `TraceSourceType == "ExponentialDistribution"`. |

### 5.12 `PeaksUsedForActivityFitting.Peaks[]`

This is the **inverted view** of §5.7: one entry per fitted peak, with the list of
sources contributing gammas to that peak nested inside, instead of one entry per source.
Each peak object carries all fields from §5.7 (built by `add_basic_peak_info()`), plus a
`Sources[]` array.

```json
"PeaksUsedForActivityFitting": {
  "Peaks": [
    {
      "Energy": "661.66", "Energy_keV": 661.66, ...   // §5.7 fields
      "Sources": [
        {
          "Energy": "661.66", "Energy_keV": 661.66,
          "HasDecayCorrection": false,
          "BranchingRatio": 0.851,
          "BranchingRatioStr": "0.851",
          "SourcePhotonsCps": 7.5e4, "SourcePhotonsCpsStr": "7.5E4",
          "SourcePhotons":   2.25e7, "SourcePhotonsStr":   "2.25E7",
          "PredictedCounts": 1.21e5, "PredictedCountsStr": "1.21E5",
          ...                                          // §5.6 source-detail fields
        }
      ]
    }
  ]
}
```

| Field (inside `Sources[]`) | Type | Units | Notes |
|---|---|---|---|
| `Energy`                | string | keV   | Gamma energy assigned to the peak from this source, formatted `%.2f`. |
| `Energy_keV`            | number | keV   | |
| `HasDecayCorrection`    | bool   |       | True if a decay-during-measurement correction was applied. |
| `SourcePhotonsCps`      | number | γ/s   | |
| `SourcePhotonsCpsStr`   | string |       | |
| `BranchingRatio`        | number | 0–1   | |
| `BranchingRatioStr`     | string |       | |
| `SourcePhotons`         | number | γ     | |
| `SourcePhotonsStr`      | string |       | |
| `PredictedCounts`       | number | cnt   | This source's contribution to the peak's predicted area. |
| `PredictedCountsStr`    | string |       | |

Plus every `add_basic_src_details()` field from §5.6.

**Note on `DecayCorrection` / `DecayCorrectionStr` in this view.** In §5.12's inverted
view, the per-source `DecayCorrection` and `DecayCorrectionStr` values are written to
the **outer peak object** rather than to each `Sources[]` element (see
`src/BatchInfoLog.cpp:1206-1207`). That means `peak.DecayCorrection` /
`peak.DecayCorrectionStr` reflects whichever source was processed *last* and is not
reliably associated with any particular source in `peak.Sources[]`. Templates that need
per-source decay corrections should iterate `Sources[].PeaksThisNucContributesTo[].ThisNucsGammasForPeak[]`
(§5.8) instead, where the correction is correctly placed on each gamma.

### 5.13 `PeakToModelComparison.UsedPeaks[]`

Short per-peak goodness-of-fit summary, suitable for a one-line-per-peak table.

| Field | Type | Notes |
|---|---|---|
| `Energy`                     | string | `printFixed(..., 2)` of the peak energy in keV. |
| `NumSigmaOff`                | string | Signed σ-distance between observed and predicted areas. |
| `ObservedOverExpected`       | string | Ratio. |
| `ObservedOverExpectedUncert` | string | Propagated uncertainty on the ratio. |

### 5.14 Chart-asset fields

These are minified JS / CSS source strings that a template embeds directly into
`<script>` / `<style>` tags to produce an interactive D3 chart. They are present
whenever InterSpec is built with D3 chart support (which is the standard configuration).

| Field | Type | Notes |
|---|---|---|
| `D3_JS`                    | string | D3 v3 minified. Embed once per page. |
| `SpectrumChart_JS`         | string | `SpectrumChartD3` implementation. |
| `SpectrumChart_CSS`        | string | `SpectrumChartD3` styles. |
| `ShieldingSourceFitPlot_JS`  | string | `ShieldingSourceFitPlot` (model-vs-observed scatter) implementation. |
| `ShieldingSourceFitPlot_CSS` | string | `ShieldingSourceFitPlot` styles. |
| `ShieldingSourceFitPlotData` | object | Pre-built data payload for `ShieldingSourceFitPlot`. Only present in act/shield reports, not in peak-fit. |

See §10.8 for a wiring snippet.

### 5.15 Application metadata

| Field | Type | Notes |
|---|---|---|
| `InterSpecCompileDate`     | string | `__DATE__` style (e.g. `"May 20 2026"`). |
| `InterSpecCompileDateIso`  | string | `YYYYMMDD` integer rendered as a string (e.g. `"20260520"`) — **not** a hyphenated ISO-8601 date. |
| `AnalysisTime`             | string | ISO-8601 *basic* format with microseconds (e.g. `"20260520T143045.123456"` — no hyphens or colons). |
| `CurrentWorkingDirectory`  | string | Working directory when the analysis ran. |
| `InterSpecExecutablePath`  | string | Absolute path to the InterSpec executable. Not present in web-deployment / iOS / Android builds. |

### 5.16 Batch-mode wrapper

When InterSpec runs the activity/shielding fit on N input files via `--batch-act-fit`,
each analysis' per-spectrum payload above appears as element `Files[i]` of a top-level
batch-summary JSON, which also carries aggregate fields.

Note that `Files` holds one entry per *analysis performed*, which is normally one per input file,
but will be more when `--multi-sample-handling each` splits a file into per-sample analyses (see
`NumInputFiles` and `NumAnalyses`).  Do not assume `Files` and `InputFiles` are the same length.

```json
{
  "Files": [
    { "Filename": "file1.n42", "Sources": [...], "Shieldings": {...}, ... },
    { "Filename": "file2.n42", "Sources": [...], "Shieldings": {...}, ... }
  ],
  "InputFiles": ["file1.n42", "file2.n42"],
  "NumInputFiles": 2,
  "NumAnalyses": 2,
  "AnalysisTime": "...",
  "CurrentWorkingDirectory": "...",
  "InterSpecCompileDate": "...",
  "InterSpecCompileDateIso": "...",
  "InterSpecExecutablePath": "...",
  "D3_JS": "...", "SpectrumChart_JS": "...", "SpectrumChart_CSS": "...",
  "ShieldingSourceFitPlot_JS": "...", "ShieldingSourceFitPlot_CSS": "..."
}
```

A template that should work in both modes can dispatch with
`{% if exists("Files") %} … {% else %} … {% endif %}` (the top-level object is
implicit, so `exists("...")` checks for a key on the root) and refer to per-file
fields as `file.Filename`, `at(file.Sources, i).Activity_uCi`, etc.

### 5.17 `SupplementalPeakInfo`

Per-peak information that supplements the fit itself: a detection-limit check of every peak, and,
for peaks that were fit but **not** used to determine the activities and shielding, the activity
each would imply on its own.

Both come from evaluating the fitted model held at its nominal answer, using the same forward model
that produced the fit — so interference between sources, self-attenuating and trace source geometry,
and decay during the measurement are all accounted for.

This object is produced by the activity/shielding computation itself, so it appears in **both** the
GUI's calculation log and the batch reports.  It is **absent entirely** when there is nothing to
report — when the fit did not converge, or when the computation is turned off (batch:
`--not-fit-peak-mda none`).  Guard with `exists("SupplementalPeakInfo")` at the top level, or
`existsIn(file,"SupplementalPeakInfo")` inside a batch-summary `Files` loop.

Two questions are answered, and they are worth keeping apart:

- **"Is there a detectable signal in this peak's region?"** — a Currie-style (ISO 11929
  gross-counts) check, made for *every* peak.  For a peak that was fit this is a quality check: a
  peak can pass the peak-fit significance tests and still sit below the level at which a signal can
  be reliably claimed.  Those show a `ShortDescription` of `"Fit, but less than Lc"`.
- **"What activity would this peak imply on its own?"** — only for peaks that were fit and were not
  used in the fit.  The useful reading is whether the unused peaks agree with the answer the used
  peaks gave; a systematic offset across all of them means something different (efficiency,
  shielding, branching ratio) than scatter does.

| Field | Type | Notes |
|---|---|---|
| `HasCurrieChecks`      | bool | Always present.  True if at least one peak's detection-limit check was computed. |
| `HasNominalActivities` | bool | Always present.  True if at least one peak has an implied activity. |
| `AnyPeakNotUsedForFit` | bool | Always present.  True if at least one fit peak was not used by the model. |
| `AllPeaks`             | object | Only present when the object itself is.  Has a `Peaks[]` array with **every** peak given to the fit. |
| `PeaksNotUsedForFit`   | object | Only present when the object itself is.  Has a `Peaks[]` array holding the subset of `AllPeaks.Peaks[]` that were fit but not used by the model.  Empty array when `AnyPeakNotUsedForFit` is false. |

Both `Peaks[]` arrays hold the same element shape (§5.17.1); the entries in
`PeaksNotUsedForFit.Peaks[]` are copies of the corresponding `AllPeaks.Peaks[]` entries, provided
separately so a report can show "peaks that disagree with the answer" without filtering in the
template.

#### 5.17.1 `SupplementalPeakInfo.*.Peaks[]` element

Peak identification — always present:

| Field | Type | Notes |
|---|---|---|
| `PeakMean`, `PeakMeanUncert`  | double | Fit peak mean and its uncertainty, keV. |
| `PeakMeanStr`                 | string | `PeakMean` pre-formatted to 5 significant figures. |
| `PeakFwhm`                    | double | FWHM in keV; for a "data defined" peak, the ROI width. |
| `PeakAmplitude`, `PeakAmplitudeUncert` | double | Observed peak area (counts) and its uncertainty. |
| `PeakAmplitudeStr`, `PeakAmplitudeUncertStr` | string | The two above, pre-formatted to 4 significant figures. |
| `HasSourceAssigned`           | bool | Whether a specific source gamma is assigned to the peak. |
| `SourceType`                  | string | `"Nuclide"`, `"X-Ray"`, `"Reaction"`, or `""`. |
| `SourceName`                  | string | E.g. `"Cs137"`; empty string if no source is assigned. |
| `SourceEnergy`                | double | Energy of the assigned source gamma, keV; `0.0` when `HasSourceAssigned` is false. |
| `SourceGammaParent`, `SourceGammaChild` | string | Decay parent/child symbols.  Only present for a nuclide source whose nuclear transition is known. |
| `UsedForFit`                  | bool | Whether this peak was used to determine the activities and shielding. |
| `IsSyntheticPeak`             | bool | True for a peak that stands in for one that was never observed — batch adds the exemplar peaks it could not fit, so the model will predict the counts each would have received.  Their `PeakAmplitude` is meaningless, and they never get an implied activity. |

Detection-limit check — always present:

| Field | Type | Notes |
|---|---|---|
| `HasCurrieCheck` | bool | Whether `CurrieCheck` holds a completed calculation.  Mirrors `CurrieCheck.CurrieComputed`. |
| `CurrieCheck`    | object | Always present.  Same shape as `Peaks[].Mda`'s counts fields — see §6.6.1 for the full key list and their presence rules. |

Model evaluation — always present:

| Field | Type | Notes |
|---|---|---|
| `ModelEvaluated`      | bool | Whether the fitted model was evaluated for this peak.  False when the peak has no nuclide assigned, its nuclide was not one of the fitted sources, it has no specific gamma assigned, or the evaluation failed. |
| `NotEvaluatedReason`  | string | Only present when `ModelEvaluated` is false.  A sentence fragment, e.g. `"Cs137 was not one of the fitted sources"`. |

Present **only when `ModelEvaluated` is true**:

| Field | Type | Notes |
|---|---|---|
| `ExpectedCounts`             | double | Counts the model predicts for this peak, from all fitted sources together. |
| `ExpectedCountsStr`          | string | The above, pre-formatted to 4 significant figures. |
| `NuclideExpectedCounts`      | double | The part of `ExpectedCounts` this peak's own nuclide contributes. |
| `OtherSourcesExpectedCounts` | double | The part all other fitted sources contribute.  This is removed before deriving an implied activity — the interference correction. |
| `SharedWithOtherSources`     | bool | Whether another fitted source contributes to this peak at all. |
| `CountsPerBq`                | double | Counts this peak receives per becquerel of its nuclide's activity.  Accounts for branching ratio at the fitted age, live time, shielding and air attenuation, detector efficiency, source geometry, and decay during the measurement. |
| `DetectorEff`                | double | Detector efficiency at this energy.  For a volumetric source, an effective value averaged over the source. |
| `ShieldingTransmission`      | double | Fraction of gammas transmitted through all shieldings. |
| `AirTransmission`            | double | Fraction transmitted through the air between the shielding and the detector; `1.0` when air attenuation is off or the geometry is fixed. |
| `IsVolumetricSource`         | bool | Whether this peak's nuclide is a self-attenuating or trace source. |

Descriptions — always present:

| Field | Type | Notes |
|---|---|---|
| `ShortDescription` | string | Table-cell sized; e.g. `"1.05x fitted"`.  Empty when no implied activity was derived. |
| `ResultSummary`    | string | A sentence with the numbers in it.  Empty when no implied activity was derived. |
| `Caveats`          | string | Notes about anything that makes the comparison less reliable; may be empty. |
| `HasCaveats`       | bool | Whether `Caveats` is non-empty.  Render caveats in their own `.MdaCaveat` span rather than burying them in the paragraph. |
| `Description`      | string | `ResultSummary` followed by `Caveats`, joined with two spaces — the caveats stay last.  Use this when you want one ready-to-print paragraph. |

Implied activity — present **only when `HasNominalActivity` is true**, which requires the peak to
have been fit, not used in the fit, and its nuclide to have a non-zero fitted activity:

| Field | Type | Notes |
|---|---|---|
| `HasNominalActivity`         | bool | Always present.  Guard every field below with it. |
| `ImpliedActivity`            | string | The activity this peak alone implies, formatted with units and any fixed-geometry postfix (e.g. `"1.23 uCi"`). |
| `ImpliedActivityUncert`      | string | Its uncertainty, likewise formatted. |
| `FitActivity`                | string | The fitted activity of this peak's nuclide, for comparison. |
| `FitActivityUncert`          | string | Its uncertainty. |
| `<any of the above>_bq`, `_kBq`, `_MBq`, `_ci`, `_mCi`, `_uCi` | double | The same activity as a plain number in that unit — e.g. `ImpliedActivity_bq`.  Use these for arithmetic and for CSV columns. |
| `ActivityPostFix`            | string | Postfix carried by the formatted strings for fixed-geometry detector responses, e.g. `"/cm2"`.  Empty in the usual case. |
| `RatioToFitActivity`         | double | Observed counts divided by the counts the model predicts.  1.0 means perfect agreement. |
| `RatioToFitActivityUncert`   | double | Its uncertainty. |
| `RatioToFitActivityStr`      | string | `RatioToFitActivity` pre-formatted to 4 significant figures. |
| `NumSigmaOff`                | double | How many standard deviations the observed counts sit from the model prediction. |
| `NumSigmaOffStr`             | string | The above, pre-formatted to 3 significant figures. |

Note that an activity family member is **omitted entirely** when the value works out non-positive —
`printToBestActivityUnits` renders those as an unreadable string of femtocuries.  Guard individual
members with `existsIn(peak,"ImpliedActivityUncert")` when a zero uncertainty is possible.

### 5.18 `PullTrend`

An automatic diagnosis of the per-peak pulls ("chi" values, i.e. `(observed-model)/uncertainty`)
as a function of energy.  The physics: if the model has the wrong *amount* of shielding, the
attenuation-vs-energy mismatch tilts the pulls monotonically with energy; if it has the wrong
effective *atomic number*, the photoelectric term (~Z⁴·⁵/E³) versus the nearly Z-independent
Compton term bends them (a low-energy "hook") even after the fit compensates the overall amount.
A low-order polynomial is fit to the pulls and the coefficient t-statistics are compared against
calibrated thresholds.  The GUI shows `Message` under the χ² chart; the bundled templates put it
in the same place.

```json
"PullTrend": {
  "HasConclusion": true,
  "Message": "Peak residuals probably indicate too much modeled shielding.",
  "SlopeT": 3.72,
  "CurvatureT": 0.41,
  "ReducedChi2": 1.18,
  "NumPeaksUsed": 12,
  "AtomicNumberDiscriminationPower": 4.83,
  "AtomicNumberDiscriminable": true
}
```

| Field | Type | Notes |
|---|---|---|
| `HasConclusion`  | bool   | A conclusion was reached *and* is worth reporting. **Always emitted.** False both when the analysis could not run at all and when it ran but found nothing conclusive — guard `Message` with this. |
| `Message`        | string | Plain-English one-sentence conclusion, e.g. `"Peak residuals probably indicate too much modeled shielding."`, `"Peak residuals may indicate the shielding's effective atomic number is too low."`, or `"A residual energy trend remains despite fitting shielding; the distance, geometry, source age, or detector response may be slightly off."`  **Always emitted**; empty string when `HasConclusion` is false.  Not localized (the GUI renders its own localized wording from the same result). |
| `SlopeT`         | number | t-statistic of the fitted linear term.  Large magnitude ⇒ a monotonic pull trend with energy ⇒ the modeled shielding *amount* is likely wrong.  Sign indicates direction. |
| `CurvatureT`     | number | t-statistic of the fitted quadratic term.  Large magnitude ⇒ the pull *shape* is wrong ⇒ the effective atomic number is likely wrong.  Only asserted as a Z conclusion when `AtomicNumberDiscriminable` is true. |
| `ReducedChi2`    | number | Reduced χ² of the trend fit itself (a measure of how well the low-order trend describes the pulls).  Values > 1 mean scatter beyond the trend; the coefficient uncertainties are already inflated by `sqrt(max(1,ReducedChi2))` so the t-statistics stay calibrated. |
| `NumPeaksUsed`   | int    | Number of pull points that survived filtering and entered the trend fit. |
| `AtomicNumberDiscriminationPower` | number | How well **this dataset** *could* resolve a reference-sized (1.5×) atomic-number error — a property of the peak energies, significances and model shielding, not of the observed pulls.  Above ~250 keV attenuation is Compton-dominated and nearly Z-independent, so only low-energy, high-significance peaks carry Z information.  0 when not computed. |
| `AtomicNumberDiscriminable`       | bool   | `AtomicNumberDiscriminationPower` cleared the calibrated gate (~2.5).  **When false, a null Z result means "could not tell", not "atomic number is fine"** — and any curvature that does appear is suppressed rather than reported as a Z signature.  Worth saying so explicitly in a report. |

> **Guarding.**  `HasConclusion` and `Message` are always present.  **The other six fields are
> present only when the trend analysis actually ran** (it is skipped, and may fail, e.g. with too
> few peaks).  Use `existsIn(PullTrend,"SlopeT")` before touching any of them:
>
> ```
> {% if exists("PullTrend") and PullTrend.HasConclusion %}{{ PullTrend.Message }}{% endif %}
> {% if exists("PullTrend") and existsIn(PullTrend,"SlopeT") %}
>   (slope t={{ printFixed(PullTrend.SlopeT,1) }}, curvature t={{ printFixed(PullTrend.CurvatureT,1) }},
>    {{ PullTrend.NumPeaksUsed }} points{% if PullTrend.AtomicNumberDiscriminable %}{% else %};
>    these peaks cannot resolve the effective atomic number{% endif %})
> {% endif %}
> ```
>
> The empty-then/`{% else %}` shape above is deliberate: a `not` operator is **not** part of the
> supported subset documented in §7 — negate with an empty `{% if %}` arm instead.

### 5.19 `DetectorEff` (and its older alias `VolumetricEff`)

Which detector-efficiency model the fit actually used, for the point sources and for the volumetric
(self-attenuating / trace) sources, and why.  One resolution serves both source kinds, so the two
names always describe the same response.  Always emitted.

```json
"DetectorEff": {
  "PointModel": "detector response at the source position (near-field & off-axis correct)",
  "VolumetricModel": "EFFTRAN transfer (near-field & off-axis correct)",
  "Note": "Auto -> EFFTRAN transfer (the one attached to the detector response)",
  "HasNote": true
}
```

| Field | Type | Notes |
|---|---|---|
| `PointModel` | string | How the point sources' efficiency was evaluated: `"detector response at the source position (near-field & off-axis correct)"` (the resolved response, through a ray fan traced once per fit), `"flat-disk (solid angle x intrinsic efficiency)"`, `"fixed-geometry intrinsic efficiency"`, or `"none (solid angle of a 1 cm disk)"`. |
| `VolumetricModel` | string | The `Method` of `VolumetricEff` below; only meaningful when `HasVolumetricSource` is true. |
| `Note` / `HasNote` | string / bool | As for `VolumetricEff`. |

`VolumetricEff` is the same information under the name the reports used before 2026-09 (kept so
older templates keep rendering):

```json
"VolumetricEff": {
  "Method": "EFFTRAN transfer (near-field & off-axis correct)",
  "Note": "Auto -> EFFTRAN transfer",
  "HasNote": true
}
```

| Field | Type | Notes |
|---|---|---|
| `Method`   | string | Human-readable name of the *resolved* method — the `Auto` setting, or a method the DRF cannot honor, has already been reduced to exactly one of: `"Monte-Carlo transfer (near-field & off-axis correct)"`, `"EFFTRAN transfer (near-field & off-axis correct)"`, `"Flat-disk (solid angle x intrinsic efficiency)"`. |
| `Note`     | string | Why it resolved that way — e.g. `"Auto -> MC transfer"`, `"EFFTRAN transfer unavailable (<reason>)"`, or `"requested near-field method, but the DRF has no CeeLo response; using flat-disk"`.  Empty when the requested method was used as-is. |
| `HasNote`  | bool   | `Note` is non-empty.  Always emitted. |

> Worth surfacing in any report covering a volumetric source: flat-disk is the legacy far-field
> approximation, so a silent fallback to it can bias a near-field (source close to the detector)
> result.  `HasNote` being true is exactly the "something was not what you asked for" signal.

## 6. JSON data model — batch peak-fit

This section describes the JSON object built by
`BatchInfoLog::add_peak_fit_results_to_json()` for the batch peak-fit tool. It shares
the §5.4 / §5.5 / §5.14 / §5.15 substructures verbatim — only the unique top-level
fields are detailed here.

### 6.1 Top-level keys

```json
{
  "Filepath": "/data/spectra/sample_001.n42",
  "Filename": "sample_001.n42",
  "ParentDir": "/data/spectra",
  "Success": true,
  "HasWarnings": false,
  "Warnings": [],
  "HasSpectrum": true,
  "EnergyCalIsRefit": true,
  "EnergyCalRefitResult": {
    "Original": { ... },        // §5.5 EnergyCal shape
    "Refit":    { ... }         // §5.5 EnergyCal shape
  },
  "FitAnyPeak":      true,
  "FitAllPeaks":     true,
  "AnyNotFitPeakMda": false,
  "ExemplarHasPeaks": true,
  "foreground":  { ... },       // §5.4
  "background":  { ... },       // §5.4 (only if a background was loaded)
  "FitPeaks":     { ... },      // §6.4
  "NotFitPeaks":  { ... },      // §6.4 (peaks the exemplar had that did not converge)
  "ExemplarPeaks":{ ... },      // §6.4 (the exemplar's original peak set)
  "PeakFitOptions": { ... },    // §6.2
  "ExemplarFile":          "/data/exemplars/cs137_template.n42",
  "ExemplarSampleNumbers": [1],
  "D3_JS": "...", "SpectrumChart_JS": "...", "SpectrumChart_CSS": "...",
  "InterSpecCompileDate":    "...",  "InterSpecCompileDateIso":  "...",
  "AnalysisTime":            "...",  "CurrentWorkingDirectory":  "...",
  "InterSpecExecutablePath": "..."
}
```

| Field | Type | Notes |
|---|---|---|
| `Filepath`              | string  | Full path to the input spectrum. |
| `Filename`              | string  | Leaf name identifying this analysis, and matching the output files written for it. Equals the input file's leaf name, unless multi-sample handling split or summed the file, in which case it carries the same `_sampleN` / `_summed` infix the output files do. |
| `SourceFilename`        | string  | The input file's leaf name, never decorated. |
| `AnalysisLabel`         | string  | Human-readable label for this analysis, e.g. `foo.n42 (sample 3)`. |
| `IsSplitFromMultiSampleFile` | bool | True when this analysis is one of several taken from a single input file. |
| `ForegroundSampleNumbers` | array[int] | Sample numbers used as the foreground. Only present when they were explicitly determined; more than one means they were summed. |
| `ParentDir`             | string  | Directory. |
| `Success`               | bool    | Fit succeeded for this file. |
| `HasWarnings`           | bool    | |
| `Warnings`              | array[string] | |
| `HasSpectrum`           | bool    | The foreground spectrum loaded successfully. |
| `EnergyCalIsRefit`      | bool    | The energy calibration was refit during peak fitting. |
| `EnergyCalRefitResult.Original` | object | The original `EnergyCal` (§5.5). Only present if `EnergyCalIsRefit`. |
| `EnergyCalRefitResult.Refit`    | object | The refit `EnergyCal` (§5.5). Only present if `EnergyCalIsRefit`. |
| `FitAnyPeak`            | bool    | At least one exemplar peak converged in this file. |
| `FitAllPeaks`           | bool    | Every exemplar peak converged. When false, `NotFitPeaks` is present. |
| `AnyNotFitPeakMda`      | bool    | Detection limits were computed for at least one peak that did not converge. |
| `ExemplarHasPeaks`      | bool    | The exemplar had any peaks at all. |
| `ExemplarFile`          | string  | Path to the exemplar file. |
| `ExemplarSampleNumbers` | array[int] | Sample numbers used from the exemplar. |
| `foreground` / `background` | object | See §5.4. |

### 6.2 `PeakFitOptions`

| Field | Type | Notes |
|---|---|---|
| `FitAllPeaks`              | bool    | Use exemplar peaks AND fit any additional peaks found. |
| `RefitEnergyCal`           | bool    | Refit the energy calibration during peak fitting. |
| `UseExemplarEnergyCal`     | bool    | Adopt the exemplar's energy calibration before fitting. |
| `WriteN42WithResults`      | bool    | Save an N42 file with the fitted peaks. |
| `ShowNonFitPeaks`          | bool    | Include exemplar peaks that did not converge in reports. |
| `OutputDir`                | string  | Output directory for generated files. |
| `CreateCsvOutput`          | bool    | Also write a per-file CSV. |
| `CreateJsonOutput`         | bool    | Also write a per-file JSON dump. |
| `OverwriteOutputFiles`     | bool    | Allow overwriting existing outputs. |
| `MultiSampleHandling`      | string  | `Auto`, `EachSampleSeparately`, or `SumAllSamples`; how input files holding more than one candidate foreground spectrum were treated. |
| `BackgroundSubFile`        | string  | Path to a background spectrum. Only present if background subtraction was requested. |
| `BackgroundSubSamples`     | array[int] | Sample numbers from the background. Only present if used. |
| `UsedExistingBackgroundPeak` | bool  | The background already had peaks that were reused. Only present if used. |
| `UseExemplarEnergyCalForBackground` | bool | Apply the exemplar's energy cal to the background. Only present if used. |
| `ReportTemplateIncludeDir` | string  | Custom template include directory (from `--report-template-include-dir`). |
| `ReportTemplates`          | array[string] | Per-file template names requested on the CLI. |
| `ReportSummaryTemplates`   | array[string] | Summary template names requested on the CLI. |
| `PeakStatThreshold`        | number  | Statistical-significance threshold for accepting peaks. |
| `PeakShapeThreshold`       | number  | Shape (χ²) threshold for accepting peaks. |
| `NotFitPeakMdaMethod`      | string  | `none`, `currie`, or `currie-and-decon`; how detection limits were computed for exemplar peaks that did not converge (§6.6.1). |
| `MdaConfidenceLevel`       | number  | Confidence level used for those limits, as a fraction; e.g. `0.95`. |
| `MdaConfidenceLevelPercent`| number  | The same, in percent; e.g. `95.0`. |
| `MdaConfidenceLevelStr`    | string  | Display form; e.g. `"95%"`, or `"1-5.7E-07"` for very high levels. |
| `MdaNumSideChannels`       | int     | Channels on each side of the peak region used to estimate the continuum. |
| `MdaRoiNumFwhm`            | number  | Width of the peak region, in multiples of the peak FWHM. **Only present in builds with `ALLOW_SPECIFY_MDA_ROI_WIDTH` enabled**, which is off by default; otherwise the width is fixed at 2.5 (the ISO 11929:2010 recommendation of mean ±1.25 FWHM, and what the GUI detection-limit tools use). The width actually used is always available per-peak as `Mda.RoiWidth_keV`. |
| `MdaSignalFractionInRoi`   | number  | Fraction of a Gaussian peak inside a region that wide. Same presence rule as `MdaRoiNumFwhm`; the per-peak `Mda.SignalFractionInRoi` is always present. |

### 6.3 Batch summary wrapper

In summary mode, the per-file shape above is wrapped.  `Files` holds one entry per *analysis
performed* - normally one per input file, but more when `--multi-sample-handling each` splits a
file into per-sample analyses.  `NumInputFiles` and `NumAnalyses` give both counts; do not assume
`Files` and `InputFiles` are the same length.

```json
{
  "Files": [
    { "Filename": "f1.n42", "FitPeaks": {...}, "foreground": {...}, ... },
    { "Filename": "f2.n42", "FitPeaks": {...}, "foreground": {...}, ... }
  ],
  "InputFiles":            ["f1.n42", "f2.n42"],
  "NumInputFiles":         2,
  "NumAnalyses":           2,
  "ExemplarFile":          "exemplar.n42",
  "ExemplarSampleNumbers": [1],
  "PeakFitOptions":        { ... },
  "D3_JS": "...", "SpectrumChart_JS": "...", "SpectrumChart_CSS": "...",
  "InterSpecCompileDate": "...", "AnalysisTime": "...", ...
}
```

### 6.4 Peak collection objects (`FitPeaks`, `NotFitPeaks`, `ExemplarPeaks`)

All three have the same shape: pre-computed sort-indices, a list of continua, and a list
of peaks. The sort indices let templates iterate peaks in different orders without
needing a custom comparator (Inja does have a built-in `sort` (see §7), but it sorts on
natural ordering only — the pre-computed indices give you multi-criterion / domain-aware
orderings like "by source-energy" or "by distance-from-source-energy").

```json
"FitPeaks": {
  "PeakSortIndex_Energy_Ascend":             [3, 0, 1, 2],
  "PeakSortIndex_Energy_Descend":            [2, 1, 0, 3],
  "PeakSortIndex_Isotope_Ascend":            [...],
  "PeakSortIndex_Isotope_Descend":           [...],
  "PeakSortIndex_Mean_Ascend":               [...],
  "PeakSortIndex_Mean_Descend":              [...],
  "PeakSortIndex_Amp_Ascend":                [...],
  "PeakSortIndex_Amp_Descend":               [...],
  "PeakSortIndex_Fwhm_Ascend":               [...],
  "PeakSortIndex_Fwhm_Descend":              [...],
  "PeakSortIndex_SrcEnergy_Ascend":          [...],
  "PeakSortIndex_SrcEnergy_Descend":         [...],
  "PeakSortIndex_RoiCounts_Ascend":          [...],
  "PeakSortIndex_RoiCounts_Descend":         [...],
  "PeakSortIndex_DistSrcEnergyToMean_Ascend":  [...],
  "PeakSortIndex_DistSrcEnergyToMean_Descend": [...],
  "PeakSortIndex_UseForActivity_Ascend":     [...],
  "PeakSortIndex_UseForActivity_Descend":    [...],
  "PeakSortIndex_UseForEnergyCal_Ascend":    [...],
  "PeakSortIndex_UseForEnergyCal_Descend":   [...],
  "Continua": [ ... ],  // §6.5
  "Peaks":    [ ... ]   // §6.6
}
```

Each `PeakSortIndex_*` is an array of zero-based indices into `Peaks[]`, sorted by the
named criterion in the named direction.

`NotFitPeaks` additionally has `"HasMdas"` (bool) — whether detection limits were computed
for the peaks it holds; see §6.6.1. `FitPeaks` and `ExemplarPeaks` never have it.

### 6.5 `Continua[]` element

Each entry describes one peak continuum (one ROI). Multiple peaks may share a single
continuum; `PeakIndexes[]` lists which `Peaks[]` entries do so.

| Field | Type | Units | Notes |
|---|---|---|---|
| `PeakIndexes`             | array[int]      |     | Indices into `Peaks[]` that use this continuum. |
| `ContinuumIndex`          | int             |     | Self-index (= the position in `Continua[]`). |
| `ContinuumType`           | string          |     | `"NoOffset"`, `"Constant"`, `"Linear"`, `"Quardratic"` (yes, this exact spelling — preserved for backwards compatibility), `"Cubic"`, `"FlatStep"`, `"LinearStep"`, `"BiLinearStep"`, `"FlatStepCDF"`, `"LinearStepCDF"`, `"BiLinearStepCDF"`, or `"External"`. |
| `NumberParameters`        | int             |     | Number of continuum parameters. |
| `IsStepContinuum`         | bool            |     | True for the stepped continuum variants. |
| `IsPolynomial`            | bool            |     | True for the polynomial variants. |
| `IsEnergyRangeDefined`    | bool            |     | True if `LowerEnergy` / `UpperEnergy` are meaningful. |
| `LowerEnergy`             | number          | keV | ROI lower edge. |
| `UpperEnergy`             | number          | keV | ROI upper edge. |
| `HasChannelRange`         | bool            |     | True if channel-level extents are available. |
| `LowerChannel`            | number          |     | Lower channel as a float. Only if `HasChannelRange`. |
| `UpperChannel`            | number          |     | Upper channel as a float. Only if `HasChannelRange`. |
| `LowerChannelInt`         | int             |     | Lower channel as integer. Only if `HasChannelRange`. |
| `UpperChannelInt`         | int             |     | Upper channel as integer. Only if `HasChannelRange`. |
| `NumberChannels`          | number          |     | Channel count as float. Only if `HasChannelRange`. |
| `NumberChannelsInt`       | int             |     | Channel count as integer. Only if `HasChannelRange`. |
| `ChannelNumbers`          | array[int]      |     | Per-channel channel numbers. Only if `HasChannelRange`. |
| `ChannelEnergies`         | array[number]   | keV | Per-channel lower edges. Only if `HasChannelRange`. |
| `ChannelContinuumArea`    | array[number]   | cnt | Per-channel continuum integral. Only if `HasChannelRange`. |
| `ContinuumArea`           | number          | cnt | Total continuum integral over the ROI. |
| `ParameterReferenceEnergy`| number          | keV | Reference energy for the continuum-parameter polynomial. |
| `Parameters`              | array[number]   |     | Continuum coefficients (in the parameter order natural to `ContinuumType`). |
| `ParameterUncertainties`  | array[number]   |     | Paired with `Parameters`. |
| `ParameterIsForFitting`   | array[bool]     |     | Per-parameter "was this a free parameter?" |
| `Warnings`                | array[string]   |     | Continuum-level warnings. Only if non-empty. |

### 6.6 `Peaks[]` element

| Field | Type | Units | Notes |
|---|---|---|---|
| `ContinuumIndex`         | int     |     | Index into the parent's `Continua[]`. Use `at(FitPeaks.Continua, peak.ContinuumIndex)`. |
| `SkewType`               | string  |     | `"NoSkew"`, `"ExGauss"`, `"DoubleExGauss"`, `"GaussPlusExGauss"`, `"GaussExp"`, `"CrystalBall"`, `"ExpGaussExp"`, `"DoubleSidedCrystalBall"`, or `"VoigtPlusExGauss"`. |
| `NumSkewParameters`      | int     |     | |
| `PeakMean`               | number  | keV | Fitted centroid. |
| `PeakMeanUncert`         | number  | keV | |
| `PeakAmplitude`          | number  | cnt | Fitted peak area. |
| `PeakAmplitudeUncert`    | number  | cnt | |
| `DataDefined`            | bool    |     | True for data-defined peaks (rare; area = data in ROI). |
| `GaussianDefined`        | bool    |     | True for Gaussian peaks (the common case). |
| `Chi2Dof`                | number  |     | χ² per degree of freedom over the ROI. **Always emitted**, but the value is only meaningful when `HasChi2Dof=true` — branch on that flag. |
| `HasChi2Dof`             | bool    |     | True if `Chi2Dof` is meaningful for this peak. |
| `PeakSigma`              | number  | keV | Gaussian σ; 0.0 for data-defined peaks. |
| `PeakSigmaUncert`        | number  | keV | |
| `PeakFwhm`               | number  | keV | `2.3548 · PeakSigma`. |
| `PeakFwhmUncert`         | number  | keV | |
| `LowerEnergy`            | number  | keV | ROI lower edge (same as continuum). |
| `UpperEnergy`            | number  | keV | ROI upper edge. |
| `RoiWidth`               | number  | keV | `UpperEnergy − LowerEnergy`. |
| `HasChannelRange`        | bool    |     | |
| `LowerChannel`           | number  |     | Only if `HasChannelRange`. |
| `UpperChannel`           | number  |     | Only if `HasChannelRange`. |
| `LowerChannelInt`        | int     |     | Only if `HasChannelRange`. |
| `UpperChannelInt`        | int     |     | Only if `HasChannelRange`. |
| `NumberChannels`         | number  |     | Only if `HasChannelRange`. |
| `NumberChannelsInt`      | int     |     | Only if `HasChannelRange`. |
| `PeakMeanChannel`        | number  | chnl | Channel-unit equivalent of `PeakMean` (fractional channel of the centroid). Only if `HasChannelRange`. |
| `PeakMeanUncertChannel`  | number  | chnl | Channel-unit equivalent of `PeakMeanUncert`. Only if `HasChannelRange`. |
| `PeakSigmaChannel`       | number  | chnl | Channel-unit equivalent of `PeakSigma`; 0.0 for data-defined peaks. Only if `HasChannelRange`. |
| `PeakSigmaUncertChannel` | number  | chnl | Channel-unit equivalent of `PeakSigmaUncert`. Only if `HasChannelRange`. |
| `PeakFwhmChannel`        | number  | chnl | Channel-unit equivalent of `PeakFwhm`. Only if `HasChannelRange`. |
| `PeakFwhmUncertChannel`  | number  | chnl | Channel-unit equivalent of `PeakFwhmUncert`. Only if `HasChannelRange`. |
| `AreaBetweenContinuumAndData` | number | cnt | Direct integration of (data − continuum) over the ROI; an alternative-to-Gaussian peak area. |
| `UseForEnergyCal`        | bool    |     | User flagged this peak for energy-cal use. |
| `UseForActivityFit`      | bool    |     | User flagged this peak for activity-fit use. |
| `UseForIsotopicsFromPeaks` | bool  |     | User flagged this peak for isotopics-from-peaks use. |
| `UseForDetEffFit`        | bool    |     | User flagged this peak for DRF-efficiency fitting. |
| `UseForDetFwhmFit`       | bool    |     | User flagged this peak for DRF-FWHM fitting. |
| `HasPeakUserLabel`       | bool    |     | |
| `PeakUserLabel`          | string  |     | User-supplied label. |
| `HasSourceAssigned`      | bool    |     | A nuclide / x-ray / reaction was assigned. |
| `SourceEnergy`           | number  | keV | Source gamma energy (0.0 if no source). |
| `SourceType`             | string  |     | `"Nuclide"`, `"X-Ray"`, `"Reaction"`, or `""`. |
| `SourceName`             | string  |     | Nuclide symbol, element name, or reaction string. |
| `SourceGammaParent`      | string  |     | Parent nuclide in a decay-chain transition. Only present if applicable. |
| `SourceGammaChild`       | string  |     | Daughter nuclide. Only present if applicable. |
| `CoefficientValues`      | array[number] | | Raw peak coefficients (Mean, Sigma, Amplitude, …, skew parameters). |
| `CoefficientUncerts`     | array[number] | | Paired with `CoefficientValues`. |
| `CoefficientFit`         | array[bool]   | | Per-coefficient "was this fit?" |
| `CoefficientNames`       | array[string] | | Human-readable names paired with `CoefficientValues`. |
| `PeakColor`              | string  |     | CSS colour (e.g. `"rgb(255,0,0)"`). |

The `Peak*Channel` fields are the channel-number equivalents of the keV `PeakMean` /
`PeakSigma` / `PeakFwhm` (and their uncertainties), computed from the spectrum's energy
calibration. They are emitted only when `HasChannelRange` is true; the width fields are
obtained by mapping the mean-centered energy interval to a channel count, so they remain
correct under non-linear calibrations.

### 6.6.1 `Peaks[].Mda` — detection limits for peaks that were not fit

When limits were computed, peaks inside `NotFitPeaks` and `ExemplarPeaks` carry a `HasMda` bool
and — when there is a limit for that particular peak — an `Mda` object. With
`--not-fit-peak-mda none` **neither key is present at all**, so guard with
`existsIn(peak,"Mda")` rather than `peak.HasMda`, or check the collection's `HasMdas` (always
present) first. Peaks inside `FitPeaks` never have either key, so templates written before this
existed are unaffected.

`NotFitPeaks` holds only the exemplar peaks that could **not** be fit; `ExemplarPeaks` holds
**every** exemplar peak, so its limits include the peaks that were fit. Use `Mda.PeakWasFit` to
tell the two cases apart, and note that a peak that was fit never gets a deconvolution limit
(§6.6.3) — for it you already have a measured area with an uncertainty.

For a peak that was not fit, the limit answers "if this exemplar peak really isn't there, how much
*could* be there and still have escaped the fit?". For a peak that **was** fit it is a quality
check: a peak can pass the peak-fit significance tests and still sit below the level at which a
signal can be reliably claimed, which shows up as a `ShortDescription` of
`"Fit, but less than Lc"`.

The counts quantities are present whenever `CurrieComputed` is true; the activity quantities only
for activity/shielding fits whose nuclide was one of the fitted sources (§6.6.2).

```
{% if existsIn(peak,"Mda") %}{{ peak.Mda.ShortDescription }} — {{ peak.Mda.Description }}{% endif %}
```

Two ready-made English strings are provided: `ShortDescription` is a brief phrase for a table
cell (`"Less than Lc"`), and `Description` is a complete report-ready paragraph. Prefer them over
assembling your own from the numbers. Branch on `ResultType`, which is a stable token, rather
than on either string.

| Field | Type | Notes |
|---|---|---|
| `ResultType`               | string  | `NotDetected` (below the decision threshold — an upper limit only), `Detected` (signal is above the decision threshold even though no peak converged), `Deficit` (significantly fewer counts than the continuum predicts), or `Error` (the limit could not be computed). |
| `ShortDescription`         | string  | Brief phrase for the result, short enough for a table cell: `"Less than Lc"`, `"Fit, but less than Lc"`, `"Greater than Lc"`, `"Fewer counts than expected"`, `"No counts in region"`, or `"Not computed"`. |
| `Description`              | string  | Ready-to-include English summary: the result, then the activity statement when applicable, then any caveats. Equals `ResultSummary` + `ActivitySummary` + `Caveats`, joined by two spaces. |
| `ResultSummary`            | string  | Just the counts statement — the first part of `Description`. |
| `ActivitySummary`          | string  | Just the activity statement, or why there isn't one; empty for plain peak fits. |
| `Caveats`                  | string  | Just the caveats — anything making the limit less reliable. Empty when there are none. |
| `HasCaveats`               | bool    | Whether `Caveats` is non-empty. The shipped HTML templates use this to render the caveats in a `.MdaCaveat` span, set apart from the result; do the same rather than burying a warning mid-paragraph. |
| `ConfidenceLevel`          | number  | Confidence level as a fraction; e.g. `0.95`. |
| `ConfidenceLevelPercent`   | number  | The same, in percent. |
| `ConfidenceLevelStr`       | string  | Display form; e.g. `"95%"`. |
| `NumSideChannels`          | int     | Channels each side of the peak region used for the continuum estimate. |
| `SignalFractionInRoi`      | number  | Fraction of the peak that lies inside the peak region; `0.9968` at the default region width. **The limit is not corrected for the signal outside the region**, so it is optimistic by roughly this factor. Negligible at the default, but a 1.0 FWHM wide region holds only 76% of the peak. When this drops below 0.99, `Description` says so in words. |
| `CurrieComputed`           | bool    | Whether the Currie-style limit succeeded. |
| `CurrieError`              | string  | Why it failed. Only present when `CurrieComputed` is false. |
| `OverlapsFitPeak`          | bool    | A peak that *was* fit falls in the evaluated region, so the limit may be biased. |
| `OverlapsOtherNotFitPeak`  | bool    | Another not-fit exemplar peak falls in the evaluated region. |
| `PeakWasFit`               | bool    | Whether this exemplar peak was actually fit in the spectrum. Always false inside `NotFitPeaks`; either value inside `ExemplarPeaks`. A deconvolution limit is never computed when true. |
| `FitPeakMean`, `FitPeakFwhm`, `FitPeakArea`, `FitPeakAreaUncert` | double | The matched fit peak's mean (keV), FWHM (keV), area (counts), and area uncertainty. **Only present when `PeakWasFit` is true.** The limit is evaluated at this peak's mean and width, rather than the exemplar's, since re-fitting the energy calibration can shift the mean. |
| `FitPeakMeanStr`, `FitPeakAreaStr`, `FitPeakAreaUncertStr` | string | The above pre-formatted (5, 4, and 4 significant figures). **Only present when `PeakWasFit` is true.** |
| `RegionIsEmpty`            | bool    | The peak region and its side channels hold essentially no counts, so the Gaussian statistics the method uses have broken down. `ResultType` is forced to `NotDetected` (a single stray count must not read as a detection), and only `DetectionLimit_counts` is meaningful — `DecisionThreshold_counts` and `UpperLimit_counts` are both zero. |

The remaining counts fields are present only when `CurrieComputed` is true:

| Field | Type | Notes |
|---|---|---|
| `DecisionThreshold_counts` | number  | Currie's critical level L<sub>c</sub>: net counts above which a signal may be reliably called detected. |
| `DetectionLimit_counts`    | number  | Currie's detection limit L<sub>d</sub>: the true net signal that may be expected to lead to detection. This is the "MDA" in counts. |
| `SourceCounts`             | number  | ISO 11929's *primary result*: nominal net counts observed above the continuum. May be negative, and is what `ResultType` is decided from. It is **not** the value to quote as the measured signal — see `BestEstimateCounts`. |
| `BestEstimateCounts`       | number  | ISO 11929-1:2019's *best estimate* (its Formula 44): the mean of the signal distribution after the impossible negative part has been truncated away and the rest renormalized. This is what the standard says to quote, and what the GUI shows. Never negative. The built-in templates put it in a "Best estimate (counts)" column beside the existing "Net counts" (which stays, and is the primary result) — as `value ± uncert` in one cell for HTML, and as two separate columns for CSV. It equals `SourceCounts` once the signal exceeds about four times its uncertainty, and exceeds it as the signal weakens — by a few percent near the decision threshold, and by tens of percent at very low counts. Given whether or not a detection was declared, per clause 10 NOTE 2. |
| `BestEstimateCountsUncert` | number  | Standard uncertainty of `BestEstimateCounts` (Formula 45). Always the smaller of it and the primary result's own uncertainty, because truncation removes part of the distribution — so do not pair one of these with the other's uncertainty. |
| `UpperLimit_counts`        | number  | Upper bound on the true signal counts, at the confidence level. |
| `LowerLimit_counts`        | number  | Lower bound. Non-negative: the coverage interval applies the same truncation at zero. |
| `PeakRegionCounts`         | number  | Gross counts in the peak region. |
| `ContinuumCounts`          | number  | Continuum counts estimated to be under the peak region. |
| `ContinuumCountsUncert`    | number  | Uncertainty on that estimate. |
| `DecisionThresholdStr`, `DetectionLimitStr`, `UpperLimitStr`, `SourceCountsStr`, `BestEstimateCountsStr`, `BestEstimateCountsUncertStr` | string | The above, pre-formatted to 4 significant figures. |
| `RoiLowerEnergy`, `RoiUpperEnergy` | number | Peak region bounds, in keV. |
| `RoiWidth_keV`             | number  | Width of the peak region, in keV (Inja cannot subtract, so this is pre-computed). |
| `RoiLowerChannel`, `RoiUpperChannel` | int  | The same, as channel numbers. |

### 6.6.2 `Peaks[].Mda` — activity fields

`HasActivity` is always present. When it is false and an activity/shielding fit tried and
failed to convert this peak, `NoActivityReason` says why (no source nuclide assigned, the
nuclide was not one of the fitted sources, the model attributes no counts to the peak, and so on).
Plain peak-fit reports have `HasActivity` false with **no** `NoActivityReason`, since
there is no source model to convert with — so guard it with `existsIn(peak.Mda,"NoActivityReason")`.
When `HasActivity` is true, these are present:

Each activity below is the activity form of the like-named counts field in §6.6.1, so
`DetectionLimitActivity` pairs with `DetectionLimit_counts`, and so on.  The minimum detectable
activity ("MDA") is `DetectionLimitActivity`.

| Field | Type | Notes |
|---|---|---|
| `Nuclide`                  | string  | The peak's parent nuclide. |
| `DetectionLimitActivity`   | string  | Formatted L<sub>d</sub> activity — the minimum detectable activity (MDA) — including any fixed-geometry postfix. |
| `DecisionThresholdActivity`| string  | Formatted L<sub>c</sub> activity. |
| `UpperLimitActivity`       | string  | Upper bound on the activity actually present, at the confidence level. **Absent whenever it would not be positive** — a `Deficit` result, or an empty region. |
| `ObservedActivity`, `ObservedActivityLower`, `ObservedActivityUpper` | string | Only when `ResultType` is `Detected`. `ObservedActivityLower` is additionally absent when the lower limit is not positive — which, since the coverage interval is truncated at zero, now means exactly zero rather than negative. `ObservedActivity` is the activity form of the primary result; `BestEstimateActivity` is what `ActivitySummary` quotes. |
| `BestEstimateActivity`, `BestEstimateActivityUncert` | string | The activity form of `BestEstimateCounts` / `BestEstimateCountsUncert` — the ISO 11929 best estimate and its standard uncertainty. Unlike `ObservedActivity` these are given whether or not a detection was declared. `ActivitySummary` quotes the value; the uncertainty appears only in these fields, so a template wanting the ISO pair has to emit both. Both are subject to the non-positive rule below, so an empty region drops them. |
| `<any of the above>_bq`, `_kBq`, `_MBq`, `_ci`, `_mCi`, `_uCi` | number | Numeric forms of each activity above; e.g. `DetectionLimitActivity_bq`. |
| `ActivityPostFix`          | string  | Postfix for fixed-geometry detector responses; e.g. `"/cm2"`. Usually empty. |
| `GammasPerBq`              | number  | Counts expected per becquerel: `activity = counts / GammasPerBq`. |
| `ShieldingTransmission`    | number  | Fraction through the fitted shielding, in (0,1]. |
| `AirTransmission`          | number  | Fraction through the intervening air; 1.0 if not modelled. |
| `DetectorEff`              | number  | Detector efficiency at the fit distance. |
| `LiveTime_s`               | number  | Spectrum live time, in seconds. |

Any activity that works out non-positive is omitted rather than printed, so **every** activity key
in this table needs an `existsIn(peak.Mda, "…")` guard, not just a `HasActivity` check.

`GammasPerBq` comes from evaluating the fitted shielding/source model at this peak's energy — the
same forward model that produced the fit — so it accounts for interference between sources,
volumetric source geometry, and decay during the measurement.  `ShieldingTransmission`,
`AirTransmission` and `DetectorEff` are the pieces of that conversion, reported for transparency;
for a volumetric source they are effective values averaged over the source rather than a single
line of sight, and `Description` notes when that is the case.  §5.17.1's `CountsPerBq` is the same
quantity, reported per peak rather than per limit.

### 6.6.3 `Peaks[].Mda` — deconvolution fields

`DeconAttempted` is always present, and is true only when the `currie-and-decon` method was
selected. This method fits the expected peak shape at a range of source strengths, rather
than just summing counts in a box.

| Field | Type | Notes |
|---|---|---|
| `DeconComputed`        | bool    | Whether it succeeded. |
| `DeconError`           | string  | Why the limit is qualified or absent. Present whenever the calculation recorded a problem — a failed computation, a scan that did not bracket an upper limit, **or** a scan whose confidence crossing did not reproduce when it was recomputed. Not restricted to `DeconComputed` being false. |
| `DeconQuantityIsCounts`| bool    | True for plain peak fits (the limit is peak counts); false for activity/shielding fits (the limit is source activity). |
| `DeconFoundUpperLimit` | bool    | Whether an upper limit was bracketed. Present when `DeconComputed` is true **and** a result object was retained (in practice the same thing). |
| `DeconBestChi2`, `DeconBestQuantity` | number | The statistic at the best fit, and the best-fit source strength. Only when `DeconComputed` is true. |
| `DeconLimitChi2`       | number  | The statistic at the limit. Present whenever `DeconComputed` is true, but **only meaningful when `DeconFoundUpperLimit` is true** — otherwise it holds a stale or sentinel value. Guard on `DeconFoundUpperLimit` before printing it. |
| `DeconContinuumNorm`   | string  | How the continuum under the peak was determined: `Floating`, `FixedByEdges`, the deprecated `FixedByFullRange`, `mixed` if regions of interest disagree, or `""` if there were no regions of interest at all. Only when `DeconComputed` is true. |
| `DeconMeasurementModel`| string  | `CurrentSpectrum` (a bound on the spectrum that was measured) or `BackgroundReference` (a prediction for a future measurement). Only when `DeconComputed` is true. |
| `DeconLimitType`       | string  | `OneSidedUpperLimit` or `CentralInterval` — different products with different thresholds. Only when `DeconComputed` is true. |
| `DeconIsPredictedSensitivity` | bool | True when the number describes a *future* measurement rather than the one in hand. **Must not be worded as an upper limit on the measured spectrum.** Only when `DeconComputed` is true. |
| `DeconSampleExposure_s`| number  | **Real** time, in seconds, of the predicted future measurement — the dwell as a user enters and reads it. Only when `DeconIsPredictedSensitivity` is true **and** the exposure is positive, so it can be absent even inside a predicted-sensitivity branch; guard with `existsIn`, not on `DeconIsPredictedSensitivity` alone. |
| `DeconOptimizerIterations`, `DeconOptimizerRestarts` | number | Continuum-optimizer health diagnostics; iterations are normally small and restarts normally zero. Only when `DeconComputed` is true and best-fit results were retained. |
| `DeconStatisticName`   | string  | Names the statistic in the χ² fields; currently always `Cash`. Same presence rule as `DeconOptimizerIterations`. |
| `DeconWarnings`        | array of string | Notes that qualify the limit without preventing it — most often that overlapping regions of interest were combined, or that the profile crossed the threshold more than twice. Each changes how the number should be read, so render them alongside the limit. Present only when non-empty. |
| `DeconUpperLimit_counts` / `DeconUpperLimitStr` | number / string | The limit, in counts. Only when `DeconFoundUpperLimit` is true **and** `DeconQuantityIsCounts` is true. |
| `DeconUpperLimit`, `DeconUpperLimit_bq`, `_ci`, `_uCi` | string / number | The limit, as an activity. Only when `DeconFoundUpperLimit` is true **and** `DeconQuantityIsCounts` is false. |
| `DeconLimitText`, `DeconBestChi2Text` | string | Pre-formatted summaries that **contain HTML** — use them only in HTML reports. Same presence rule as `DeconUpperLimit` (activity mode). |
| `MethodsDisagree`      | bool    | `DeconOverCurrieRatio` falls outside the empirical agreement window of [0.66, 1.23]. A model-discrepancy warning (usually a curved continuum, or an interfering peak in the region), not a statistical test. Only when `DeconComputed` is true. |
| `DeconOverCurrieRatio` | number  | The deconvolution upper limit divided by the **Currie upper limit** (`UpperLimit_counts`, not L<sub>d</sub>). Present only when both are usable. |

> **The χ² fields are not chi-squares.** `DeconBestChi2` and `DeconLimitChi2` hold the **Cash
> statistic** — a Poisson deviance, `2·Σ(E − n + n·ln(n/E))`, summed once per channel (Cash 1979).
> The name is kept for continuity with the display and with older reports. Only *differences* of it
> between trial source strengths are meaningful; a difference follows the same asymptotic χ² law a
> Gaussian one would, which is what sets the confidence threshold. Its absolute value does **not**
> test goodness of fit: E[statistic] is not (channels − parameters), so a χ²/DOF built from these
> is not a fit-quality measure and should not be presented as one.

The two methods look at the same data, so a large disagreement says something about the data
rather than about the methods — most often a curved continuum under the peak, or an interfering
peak in the region, i.e. exactly where the gross-counts assumption breaks down. When
`MethodsDisagree` is true the reason is also appended to `Caveats` (and therefore `Description`),
so simply rendering those picks it up.

Every field above other than `DeconAttempted` needs a guard; the safe pattern is
`{% if peak.Mda.DeconAttempted %}{% if peak.Mda.DeconComputed %}{% if peak.Mda.DeconFoundUpperLimit %}…`.
That pattern is **not sufficient on its own** for a key carrying an extra presence condition —
notably `DeconSampleExposure_s`, which also requires a positive exposure. Referencing an absent key
is a hard `RenderError` that aborts the **entire** report, not just that sentence, so wrap it in
`{% if existsIn(peak.Mda,"DeconSampleExposure_s") %}`. The bundled templates do this.

**Times are reported in real time.** A user enters a dwell as a real time and every reported time is
a real time; the likelihood works internally in live time. A detector reporting only one of the two
times is taken to have zero dead time, making them equal.

In an activity/shielding report, peaks that could not be converted to an activity fall back to the
counts-based limit, so `DeconQuantityIsCounts` can differ from peak to peak within one report.

### 6.7 Reused common objects

The following shapes are described in §5 and used **identically** in the peak-fit
payload — refer back rather than re-document them here:

- `foreground` / `background` — §5.4.
- `EnergyCal` (also inside `EnergyCalRefitResult.{Original,Refit}`) — §5.5.
- `D3_JS` / `SpectrumChart_JS` / `SpectrumChart_CSS` — §5.14. (Note: peak-fit does **not**
  include `ShieldingSourceFitPlot_*` — that pair is act/shield-only.)
- `InterSpecCompileDate` / `AnalysisTime` / `CurrentWorkingDirectory` /
  `InterSpecExecutablePath` — §5.15.

## 7. Inja syntax — what works

Inja 3.5.0 (vendored at `external_libs/SpecUtils/3rdparty/inja/`). The features below
are confirmed to work against the bundled templates.  Two of them the detection-limit blocks now
rely on: **`{% else if %}`** (chained else-if, spelled with a space; one `{% endif %}` closes the
whole chain, not one per arm) and **string equality** (`{% if x == "SomeValue" %}` against a
literal). `set_trim_blocks(true)` is on, so
the first newline after a block is dropped.

```jinja
{# This is a comment.  Inja drops comments before parsing. #}

{{ variable }}                                    {# interpolate a top-level field #}
{{ object.field }}                                {# dotted access into nested objects #}
{{ array.0.field }}                               {# integer-literal index into an array #}
{{ at(array, n).field }}                          {# variable-index access — use this for arr[expr] #}

{% if Some.Bool %}                                {# boolean conditional #}
  ...
{% else %}
  ...
{% endif %}

{% if existsIn(obj, "OptionalKey") %} ... {% endif %}     {# guard a possibly-missing field #}
{% if exists("TopLevelKey") %} ... {% endif %}            {# guard a possibly-missing top-level key #}

{% for src in Sources %}                          {# loop over an array #}
  {{ loop.index1 }}.  {{ src.Nuclide }}: {{ src.Activity }}
{% endfor %}

{% for i in range(length(arr)) %}                 {# explicit-index loop #}
  {{ at(arr, i) }}
{% endfor %}

{% include "default-act-fit-html-results" %}      {# embed a pre-registered bundled template #}
{% include "./snippets/header.inc" %}             {# embed a file under the include dir (only when set_search_included_templates_in_files is on) #}

## for file in Files                              {# `##` line-statement syntax (equivalent to `{% for %}`) #}
  {{ file.Filename }}
## endfor

## if file.Success
  ...
## else
  ...
## endif
```

> **⚠️ Markdown footgun.** Inja treats **any** line that begins with `##` as a
> control statement — including Markdown H2/H3/H4 headings (`## Heading`,
> `### Heading`, `#### Heading`). Writing a Markdown `## Results` line in a
> `.tmplt.md` file produces `parser_error: ... expected statement, got 'Results'`.
> See §8.6 for the fix. If you are authoring a Markdown template, **prefer the
> `{% ... %}` block syntax** over `## ...` line statements and use the escape in
> §8.6 for headings.


`loop` variables inside a `for` block:

| Variable | Meaning |
|---|---|
| `loop.index`     | 0-based iteration index. |
| `loop.index1`    | 1-based iteration index. |
| `loop.is_first`  | True on first iteration. |
| `loop.is_last`   | True on last iteration. |
| `loop.parent`    | The enclosing loop's `loop` (when nesting `for` blocks). |

Built-in Inja functions usable in expressions: `length(arr)`, `range(n)` (`0`..`n-1`),
`exists("key")`, `existsIn(obj, "key")`, `upper(s)`, `lower(s)`, `default(x, fallback)`,
`at(arr, index)`, `first(arr)`, `last(arr)`, `min(a, b)`, `max(a, b)`, `round(x, digits)`,
`sort(arr)`, `join(arr, sep)`, `int(x)`, `float(x)`, `even(x)`, `odd(x)`,
`divisibleBy(x, n)`, `isString(x)` / `isArray(x)` / `isObject(x)` / `isNumber(x)` /
`isBoolean(x)` / `isFloat(x)` / `isInteger(x)`, `capitalize(s)`, `replace(s, old, new)`.

## 8. What does NOT work

This Inja version has several rough edges. Workarounds:

### 8.1 Bracket indexing with a non-literal index

```jinja
{{ Sources[loop.index].Activity }}      {# DOES NOT WORK — Inja errors out #}
{{ at(Sources, loop.index).Activity }}  {# Use this instead #}
{{ Sources.0.Activity }}                {# Integer-literal indexing IS allowed (note the `.0`, not `[0]`) #}
```

Prefer `at(arr, expr)` whenever the index isn't a literal.

### 8.2 Expression-side arithmetic is fragile

The Inja expression evaluator parses numbers via JSON and is finicky about types,
NaN / inf, and operator precedence. Things like `{{ src.Activity_uCi * 1000 }}` or
`{{ x / y }}` will sometimes succeed and sometimes throw a parse error, particularly
when an operand is null or comes from a callback. **Prefer the pre-computed unit
variants** — e.g. use `src.Activity_pCi` instead of `src.Activity_uCi * 1e6`. When you
do need a derived value, compute it in C++ and add a field, or use `printFixed` /
`printCompact` callbacks for formatting.

### 8.3 No HTML auto-escaping

Inja does not HTML-escape interpolated values. Most fields in the payload are
machine-generated and safe, but user-supplied strings (`SpectrumTitle`,
`SpectrumRemarks`, `PeakUserLabel`, `SourceName`, etc.) could in principle contain `<`
or `&`. If you embed those in an HTML report, sanitise them in C++ before adding to
the JSON, or accept that they will be rendered as raw markup.

### 8.4 Compound conditionals can surprise

Inja does support `and` / `or` / `not`, but mixing them with callbacks or array access
sometimes produces a parse error. If you hit one, split into nested `{% if %}` blocks.

### 8.5 `{% include %}` path requirements

Custom templates can `{% include "name" %}` either:

- A bundled named include (see §4 — `default-act-fit-html-results` etc.); or
- A file path **relative to the include-dir** that was passed to `get_default_inja_env()`.

When InterSpec opens a custom user template, the include-dir is set to the directory
containing the template, **with a trailing path separator**. If you author a custom
include, place it alongside the master template and reference it by relative path.

### 8.6 Markdown headings (`##`, `###`, …) collide with the `##` line-statement

Inja's default line-statement prefix is `##`, and the parser checks for it at the
start of **every** line (column 0). That means any of these Markdown lines in a
`.tmplt.md` file will fail to parse:

```markdown
## Results               ← parsed as "## <statement>", errors with: expected statement, got 'Results'
### Sub-heading          ← starts with `##`, also fails
#### Detail              ← starts with `##`, also fails
```

The error you'll see is something like:
```
parser_error: line 8, column 4 of 'my_report.tmplt.md': expected statement, got 'Fitted'.
```

You have three options:

1. **Escape the `##` via interpolation** (the cleanest fix — the rendered output is
   exactly the heading you intended):
   ```jinja
   {{ "## " }}Fitted Activities
   {{ "### " }}Sub-section
   ```
   This produces literal `## Fitted Activities` / `### Sub-section` in the output,
   which Markdown then renders as H2 / H3 normally.

2. **Use a single `#` for H1** (only `##` and longer prefixes collide, so a single
   `#` Markdown heading is safe):
   ```markdown
   # Top-level heading is fine, no escape needed
   ```

3. **Use Setext-style headings** (underlines instead of leading `#`s):
   ```markdown
   Fitted Activities
   =================

   Sub-section
   -----------
   ```
   Setext only supports H1 (`===`) and H2 (`---`); deeper levels still need option 1.

Inside `.tmplt.html` / `.tmplt.txt` / `.tmplt.csv` files, none of this is a problem
unless you intentionally write a line starting with `##` (which you almost certainly
wouldn't).

**Consequence for Markdown templates:** prefer the `{% for %} ... {% endfor %}` /
`{% if %} ... {% endif %}` block syntax over the `## for` / `## if` line-statement
syntax. Both work, but if you stick to `{% %}` you can't accidentally turn a heading
into a syntax error by adding text starting with `##`.

## 9. Custom callbacks

`BatchInfoLog::get_default_inja_env()` registers two callbacks for this tool. (The
Isotopics-by-Nuclides tool has additional callbacks like `pct`, `safe_html`,
`scientific`, etc. — those are **not** registered for ShieldSourceFitLog and must not
be used.)

| Callback | Args | Returns | Description |
|---|---|---|---|
| `printFixed(value, decimals)`     | (number, int) | string | Format `value` with fixed-point notation and exactly `decimals` digits after the decimal point. E.g. `printFixed(3.14159, 2)` → `"3.14"`. |
| `printCompact(value, sig_figs)`   | (number, int) | string | Format `value` with `sig_figs` significant figures, using whichever of fixed or scientific notation is more compact. Backed by `SpecUtils::printCompact`. E.g. `printCompact(1.23456e5, 4)` → `"1.234E5"`; `printCompact(0.0042, 3)` → `"0.00420"`. |

Prefer these callbacks over Inja expression arithmetic for any formatting that
involves a decimal point or scientific notation.

### `at(array, index)`

`at` is provided by Inja itself (not registered here), and is the recommended workaround
for bracket-indexing with a variable index (§8.1):

```jinja
{{ at(FitPeaks.Continua, peak.ContinuumIndex).LowerEnergy }}
```

## 10. Common patterns

### 10.1 Activity / Shielding — list each source's activity ± uncertainty in µCi

```jinja
## for src in Sources
  - {{ src.Nuclide }}: {{ printCompact(src.Activity_uCi, 5) }} µCi
    {% if src.ActivityIsFit %}± {{ printCompact(src.ActivityUncert_uCi, 3) }} µCi ({{ src.ActivityUncertPercent }}%){% endif %}
    {% if src.AgeIsFittable %}, age {{ printCompact(src.AgeDays, 4) }} d{% if src.AgeIsFit %} ± {{ printCompact(src.AgeUncertDays, 3) }} d{% endif %}{% endif %}
## endfor
```

### 10.2 Activity / Shielding — walk shieldings, branch on generic vs material

```jinja
## for shield in Shieldings.Shields
  Layer {{ loop.index1 }}: {{ shield.Name }}
  {% if shield.IsGeneric %}
    Generic: AN = {{ shield.AN }}, AD = {{ shield.AD }} (fit: AN={{ shield.FitAN }}, AD={{ shield.FitAD }})
  {% else %}
    Material: {{ shield.Formula }}, density {{ shield.DensityStr }}
    Thickness: {{ shield.Thickness }}{% if shield.DimensionIsFit.0 %} (fit){% endif %}
  {% endif %}
## endfor
```

### 10.3 Activity / Shielding — per-gamma counts attribution

For every nuclide, for every peak it contributes to, for every gamma line of that
nuclide contributing to that peak:

```jinja
## for src in Sources
  {{ src.Nuclide }}:
  {% for peak in src.PeaksThisNucContributesTo %}
    Peak @ {{ peak.Energy }} keV (observed {{ peak.CountsStr }} ± {{ peak.CountsUncertStr }}, predicted {{ printCompact(peak.PredictedCounts, 4) }}, σ-off {{ printCompact(peak.PredictedNumSigmaOff, 3) }})
    {% for g in peak.ThisNucsGammasForPeak %}
      - {{ g.Energy }} keV, BR={{ g.BranchingRatioStr }}, predicted {{ g.PredictedCountsStr }} cnt
    {% endfor %}
  {% endfor %}
## endfor
```

### 10.4 Peak-fit — iterate `FitPeaks.Peaks[]` in energy order

```jinja
{% for i in FitPeaks.PeakSortIndex_Energy_Ascend %}
  {% set peak = at(FitPeaks.Peaks, i) %}
  {{ printFixed(peak.PeakMean, 2) }} keV
    FWHM {{ printFixed(peak.PeakFwhm, 3) }} keV,
    area {{ printCompact(peak.PeakAmplitude, 5) }} ± {{ printCompact(peak.PeakAmplitudeUncert, 3) }},
    χ²/dof {% if peak.HasChi2Dof %}{{ printFixed(peak.Chi2Dof, 2) }}{% else %}—{% endif %},
    continuum {{ at(FitPeaks.Continua, peak.ContinuumIndex).ContinuumType }}
{% endfor %}
```

### 10.5 Peak-fit — energy-cal refit comparison

```jinja
{% if EnergyCalIsRefit %}
  Energy-cal refit ({{ EnergyCalRefitResult.Original.Type }} → {{ EnergyCalRefitResult.Refit.Type }}):
    Original coefficients: {% for c in EnergyCalRefitResult.Original.Coefficients %}{{ printCompact(c, 6) }} {% endfor %}
    Refit coefficients:    {% for c in EnergyCalRefitResult.Refit.Coefficients %}{{ printCompact(c, 6) }} {% endfor %}
{% endif %}
```

### 10.6 Batch — dispatch single-file vs `Files[]` shape

```jinja
{% if exists("Files") %}
  Batch run: {{ length(Files) }} files
  ## for file in Files
    - {{ file.Filename }}: {{ file.Success }}
  ## endfor
{% else %}
  Single file: {{ Filename }}
{% endif %}
```

### 10.7 Guarding optional fields

```jinja
{% if existsIn(foreground, "DoseRate_uSvPerHour") %}
  Ambient dose rate: {{ printFixed(foreground.DoseRate_uSvPerHour, 3) }} µSv/h
{% endif %}

{% if src.IsTraceSource and src.TraceActivityType == "ExponentialDistribution" %}
  Relaxation length: {{ src.TraceRelaxationLength_cm }} cm
{% endif %}
```

### 10.8 Embedded interactive chart

For an HTML report that shows the foreground spectrum with the shielding/source-fit
overlay:

```html
<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>{{ Filename }} — Activity/Shielding Fit</title>
  <style>{{ SpectrumChart_CSS }}</style>
  <style>{{ ShieldingSourceFitPlot_CSS }}</style>
</head>
<body>
  <h1>{{ Filename }} — χ² = {{ printFixed(FitChi2, 2) }}, dof = {{ NumDof }}</h1>

  <div id="spec-chart" style="width:100%; height:400px;"></div>
  <div id="fit-chart"  style="width:100%; height:300px;"></div>

  <script>{{ D3_JS }}</script>
  <script>{{ SpectrumChart_JS }}</script>
  <script>{{ ShieldingSourceFitPlot_JS }}</script>
  <script>
    // Spectrum chart (foreground)
    const specChart = new SpectrumChartD3("spec-chart", {});
    const specData  = {{ foreground.spectrum }};
    specChart.setData( specData );

    // Model-vs-observed fit chart
    const fitChart  = new ShieldingSourceFitPlot("fit-chart", {});
    const fitData   = {{ ShieldingSourceFitPlotData }};
    fitChart.setData( fitData );
  </script>
</body>
</html>
```

(For peak-fit reports, drop the `ShieldingSourceFitPlot*` bits — those fields are not
present in the peak-fit payload.)

## 11. Example templates

Three short, complete starting points. They are intentionally minimal; consult
`act_fit.tmplt.html`, `std_fit_log.tmplt.txt`, `peak_fit.tmplt.html`, and
`std_peak_fit_summary.tmplt.html` (shipped alongside this README) for fully-featured
production examples.

### 11.1 Activity / Shielding CSV — one row per nuclide

Filename: `act_summary.tmplt.csv`

```jinja
Filename,Nuclide,Activity_uCi,ActivityUncert_uCi,UncertPercent,Age_d,Chi2,DoF
## for src in Sources
{{ Filename }},{{ src.Nuclide }},{{ src.Activity_uCi }},{% if src.ActivityIsFit %}{{ src.ActivityUncert_uCi }},{{ src.ActivityUncertPercent }}{% else %},{% endif %},{{ src.AgeDays }},{{ printFixed(FitChi2, 3) }},{{ NumDof }}
## endfor
```

### 11.2 Activity / Shielding plain-text one-pager

Filename: `act_onepager.tmplt.txt`

```jinja
=== Activity / Shielding Fit Results ===
File:     {{ Filename }}
Detector: {{ Detector.Name }}
{% if not FixedGeometryDetector %}
Distance: {{ ActShieldFitSetup.Distance }} ({{ ActShieldFitSetup.Geometry }})
{% else %}
Geometry: fixed ({{ ActShieldFitSetup.FixedGeometryType }})
{% endif %}
Fit:      χ² = {{ printFixed(FitChi2, 2) }}, dof = {{ NumDof }}, fcn calls = {{ NumberFcnCalls }}

-- Sources --
## for src in Sources
{{ src.Nuclide }}:
  Activity: {{ src.Activity }}{% if src.ActivityIsFit %} ± {{ src.ActivityUncert }} ({{ src.ActivityUncertPercent }}%){% endif %}
  Age:      {{ src.Age }}{% if src.AgeIsFit %} ± {{ src.AgeUncert }}{% endif %}{% if not src.AgeIsFittable %} (not fittable){% endif %}
  {% if src.IsTraceSource %}Trace:    {{ src.TraceDisplayActivity }} ({{ src.TraceActivityType }}){% endif %}
  {% if src.IsSelfAttenSource %}Self-att: in layer "{{ src.SelfAttenShieldName }}", mass-frac {{ printFixed(src.SelfAttenMassFrac, 6) }}{% endif %}
## endfor

-- Shieldings ({{ Shieldings.Geometry }}, {{ Shieldings.NumberShieldings }} layers) --
## for shield in Shieldings.Shields
[{{ loop.index1 }}] {{ shield.Name }}: thickness {{ shield.Thickness }}{% if shield.DimensionIsFit.0 %} (fit){% endif %}{% if not shield.IsGeneric %}, ρ = {{ shield.DensityStr }}{% else %}, AN = {{ shield.AN }}, AD = {{ shield.AD }}{% endif %}
## endfor

-- Per-peak goodness of fit --
{% for p in PeakToModelComparison.UsedPeaks %}
  {{ p.Energy }} keV: σ-off {{ p.NumSigmaOff }}, O/E {{ p.ObservedOverExpected }} ± {{ p.ObservedOverExpectedUncert }}
{% endfor %}

Analysed at {{ AnalysisTime }} by InterSpec (compiled {{ InterSpecCompileDate }}).
```

### 11.3 Peak-fit Markdown summary

Filename: `peak_fit_summary.tmplt.md`

Note the `{{ "## " }}` and `{{ "### " }}` escapes on the heading lines — without them
Inja's `##` line-statement would try to parse the heading text as a control statement
and fail (see §8.6).

```jinja
# Peak-fit results — {{ Filename }}

- Exemplar: `{{ ExemplarFile }}` (samples {{ ExemplarSampleNumbers }})
- Spectrum: {{ foreground.LiveTime }} live ({{ foreground.NumberChannels }} channels, {{ printFixed(foreground.LowerSpectrumEnergy, 1) }} – {{ printFixed(foreground.UpperSpectrumEnergy, 1) }} keV)
- All exemplar peaks converged: **{{ FitAllPeaks }}**
{% if EnergyCalIsRefit %}- Energy calibration refit during fit.{% endif %}

{{ "## " }}Fitted peaks

| Energy (keV) | FWHM (keV) | Area | Area uncert | χ²/dof | Source |
|---:|---:|---:|---:|---:|---|
{% for i in FitPeaks.PeakSortIndex_Energy_Ascend %}
{% set peak = at(FitPeaks.Peaks, i) %}
| {{ printFixed(peak.PeakMean, 2) }} | {{ printFixed(peak.PeakFwhm, 3) }} | {{ printCompact(peak.PeakAmplitude, 5) }} | {{ printCompact(peak.PeakAmplitudeUncert, 3) }} | {% if peak.HasChi2Dof %}{{ printFixed(peak.Chi2Dof, 2) }}{% else %}—{% endif %} | {% if peak.HasSourceAssigned %}{{ peak.SourceName }} ({{ peak.SourceType }}){% else %}—{% endif %} |
{% endfor %}

{% if existsIn(NotFitPeaks, "Peaks") and length(NotFitPeaks.Peaks) > 0 %}
{{ "## " }}Peaks that did not converge
{% for p in NotFitPeaks.Peaks %}
- {{ printFixed(p.PeakMean, 2) }} keV (exemplar){% if existsIn(p, "Mda") %} — {{ p.Mda.Description }}{% endif %}
{% endfor %}
{% endif %}

_Analysed at {{ AnalysisTime }} by InterSpec ({{ InterSpecCompileDate }})._
```

## 12. Suggested LLM prompt

When asking an LLM to author a custom template using this document, use a prompt
along these lines:

> You are writing an Inja template for InterSpec's Activity/Shielding source fit
> (or batch peak-fit) tool. Use **only** fields documented in §5 (act/shield) or §6
> (peak-fit) of the README — do not invent fields. Specifically:
>
> 1. **Pick the family.** Confirm with me whether the report should consume the
>    Activity/Shielding payload (§5) or the batch peak-fit payload (§6), and whether
>    it should support both single-file mode and the batch `Files[]` wrapper.
> 2. **Pick the file extension.** `.tmplt.html` for an interactive report (you may
>    embed the D3 spectrum chart and, for act/shield, the shielding/source-fit plot
>    using §10.8). `.tmplt.txt` / `.tmplt.csv` / `.tmplt.md` for non-HTML outputs.
> 3. **Format numerics with `printFixed(value, decimals)` or `printCompact(value, sig_figs)`** —
>    do not perform expression-side arithmetic; use the pre-computed unit variants
>    instead (e.g. `Activity_uCi`, `Distance_m`).
> 4. **Guard optional fields** with `exists("...")` or `existsIn(obj, "...")` — see
>    §5/§6 "only present if …" notes.
> 5. **Use `at(arr, index)` for variable indexing** — never `arr[expr]`.
> 6. **Do NOT use** `pct()`, `safe_html()`, `scientific()` or other callbacks not
>    listed in §9 — those exist in a different InterSpec tool but are not registered
>    here.
> 7. When the report is HTML and embeds an interactive chart, base it on §10.8 and
>    embed `D3_JS`, `SpectrumChart_JS`, `SpectrumChart_CSS` (and, for act/shield,
>    `ShieldingSourceFitPlot_JS/CSS` + `ShieldingSourceFitPlotData`) once at the top.
> 8. **Markdown footgun (`.tmplt.md` only):** Inja's `##` line-statement prefix
>    collides with Markdown H2/H3/H4 headings. Any line that starts with `##` (or
>    `###`, `####`, …) is parsed as a control statement and will error out with
>    something like `expected statement, got 'Heading'`. When writing Markdown:
>    - Escape headings via interpolation: `{{ "## " }}Section` (renders as a literal
>      H2 in the output).
>    - Single `#` H1 headings are safe — no escape needed.
>    - Setext headings (`Title` on one line, `===` or `---` underneath) are also safe.
>    - For control flow inside Markdown templates, prefer `{% for %} ... {% endfor %}`
>      and `{% if %} ... {% endif %}` blocks over the `## for` / `## if` line-statement
>      forms, so a heading typo can never break parsing.
>    See §8.6 for details.
> 9. **Sanity-check the result.** Before declaring done, walk through the template
>    once: confirm every field reference appears in §5 or §6 of the README, every
>    callback is in §9, no line begins with `##` unless it is an actual control
>    statement, and any `{% include %}` target either is a §4 named include or a file
>    that will exist alongside the template.
>
> When done, tell me to save the file as `<chosen_name>.tmplt.<ext>` under
> `<writable-data-dir>/act_shield_fit_templates/` (per §3.2 for my platform — or use
> **Help → About InterSpec → Data** in the app to find and open it, per §3.5) and
> restart InterSpec to pick it up.
>
> Here is the user's request: …

## 13. Where the C++ side lives

For LLMs that have web access, the authoritative source for everything above is in:

- [`src/BatchInfoLog.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/BatchInfoLog.cpp) —
  all JSON-population functions (`shield_src_fit_results_to_json`,
  `add_peak_fit_results_to_json`, `add_basic_src_details`, `add_basic_peak_info`,
  `add_gamma_info_for_peak`, `add_hist_to_json`, `add_energy_cal_json`,
  `add_exe_info_to_json`, `add_act_shield_fit_options_to_json`, `get_default_inja_env`,
  `add_mda_to_json`, `add_currie_check_to_json`, `add_supplemental_peak_info_to_json`).
- [`InterSpec/BatchInfoLog.h`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec/BatchInfoLog.h) —
  public interfaces.
- [`src/ShieldingSourceDisplay.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/ShieldingSourceDisplay.cpp) —
  `showCalcLog()` (the GUI-side renderer; also where the writable directory is scanned).
- [`InterSpec/ShieldingSourceFitCalc.h`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec/ShieldingSourceFitCalc.h) —
  the `ModelFitResults` / `SourceFitDef` / `ShieldingInfo` / `SupplementalPeakInfo` data
  structures that the JSON payload mirrors.
- [`src/ShieldingSourceFitCalc.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/ShieldingSourceFitCalc.cpp) —
  `compute_supplemental_peak_info()`, which produces §5.17 for both the GUI and batch.
- [`InterSpec/DetectionLimitCalc.h`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec/DetectionLimitCalc.h) —
  `PeakCurrieCheck` and `currie_check_for_peak()`, the shared per-peak detection-limit check
  behind both `Mda` (§6.6.1) and `CurrieCheck` (§5.17.1).
- [`InterSpec/ShieldSourcePullTrend.h`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec/ShieldSourcePullTrend.h) —
  `TrendResult` and `conclusion_report_text()`, the pull-trend diagnosis behind §5.18 (the
  header comments cover the physics and the calibrated thresholds).
- [`src/BatchPeak.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/BatchPeak.cpp) —
  the batch peak-fit driver and `Files[]` wrapping.
- [`InterSpec/BatchPeak.h`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec/BatchPeak.h) —
  `BatchPeakFitResult` data structure.
- [Inja](https://github.com/pantor/inja) (project) and
  [Inja docs](https://pantor.github.io/inja/) — for the templating engine itself.

InterSpec uses Inja 3.5.0, vendored at `external_libs/SpecUtils/3rdparty/inja/`.
