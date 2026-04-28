# RelActAuto Report Template Authoring Guide

> **Intent.** This document is meant to be handed to an LLM together with a description of the
> report a user wants. The LLM should produce a complete Inja template (`*.tmplt.html`,
> `*.tmplt.txt`, `*.tmplt.csv`, `*.tmplt.md`, etc.) that consumes the JSON data documented below.
> Three default templates ship in `InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/`
> (`std_rel_eff_summary.tmplt.html`, `std_rel_eff_summary.tmplt.txt`, and
> `std_multi_file_summary.tmplt.html`); they are themselves examples of what the LLM should
> produce.  The InterSpec code repository can be accessed at
> https://github.com/sandialabs/InterSpec/ .
>
> **Two report flavours** are documented here:
> - **Single-file (per-spectrum) report** — sections §1 through §10 — consumes the JSON returned
>   by `RelActAutoReport::solution_to_json()` for one spectrum.  Used for both stand-alone
>   single-spectrum reports and as the per-file pages of a batch run.
> - **Multi-file batch summary report** — see §11 — consumes the JSON built by
>   `BatchRelActAuto::run_in_files()` to summarise an N-file batch (each file's `solution_to_json`
>   payload is nested under `Files[i]`).
>
> The two flavours map directly onto the InterSpec command-line and GUI as follows:
>
> | Flavour | CLI arg | GUI control (Batch dialog → Isotopics tab) | Built-in shorthand |
> |---|---|---|---|
> | **Per-file** | `--file-report-template` | "Custom Per-File" report uploader | `"html"`, `"txt"`, `"json"` |
> | **Multi-file summary** | `--summary-report-template` | "Custom Summary" report uploader | `"html-summary"` (also accepted: `"summary"`) |
>
> Pass a custom-template filename or absolute path on either flag/uploader and InterSpec will
> resolve it the same way it resolves the built-in shorthands (see §2).
>
> **Last updated:** 2026-04-27

## 1. What this is

`InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/` is the directory in
[InterSpec](https://github.com/sandialabs/InterSpec) that holds [Inja](https://github.com/pantor/inja)
templates for rendering reports of a `RelActCalcAuto::RelActAutoSolution` — the result of a
relative-efficiency / relative-activity ("isotopics by nuclides") fit. Each template consumes a
JSON payload built by `RelActAutoReport::solution_to_json()` and produces a text artifact (HTML,
plain text, CSV, Markdown, JSON, …).

The remainder of this document, sections §1 through §10, focuses on **single-file (per-spectrum)
templates** — i.e. templates that consume the JSON for one solved spectrum.  For
**multi-file batch-summary templates** (templates that consume an aggregate JSON containing all
files in a batch), see [§11 below](#11-multi-file-batch-summary-templates).

Three templates ship in `InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/` by
default:

| File | Purpose |
|---|---|
| `std_rel_eff_summary.tmplt.html` | **Per-file** full HTML report with embedded D3 spectrum chart and rel-eff chart |
| `std_rel_eff_summary.tmplt.txt` | **Per-file** plain-text summary |
| `std_multi_file_summary.tmplt.html` | **Multi-file** batch-summary HTML report (see §11) |

A user can drop their own `.tmplt.*` files into
`InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/` (or into the writable
`<InterSpec::writableDataDirectory()>/IsotopicsByNuclidesReportTmplts/` directory, or into any
directory they pass to the report API as `include_dir`) and select them at render time by passing
the filename.

### 1.1 What the "isotopics by nuclides" tool does

The tool takes a gamma-ray spectrum and, given a user-supplied list of candidate nuclides and
a user-supplied list of energy regions (ROIs) to consider, jointly fits:

- the **photo-peak amplitudes** in each ROI (Gaussian peaks on a continuum);
- one or more **relative-efficiency curves** that describe how the detector + shielding + source
  geometry attenuate gammas as a function of energy;
- the **relative activities** of each nuclide (i.e., the activity of each nuclide *expressed as
  a fraction of the total* — not absolute Bq);
- optionally, an **energy-calibration adjustment** (offset, gain, sometimes cubic) and per-pair
  **deviation pairs**, FWHM coefficients, peak-skew parameters, and per-nuclide **ages**.

A fit may fail for a bad configuration or no numerical convergence. The output is the `RelActAutoSolution` C++ struct that `solution_to_json()` serialises to JSON, that the Inja templating engine then applies to the report template that this file describes how to create.

For uranium and plutonium fits the relative activities directly give:

- **U-235 enrichment** (mass fraction of U-235 within total uranium);
- **Pu isotopic mass fractions** (Pu-238/239/240/241/242), with optional Pu-242-by-correlation
  estimation, Pu-aging, and "back-decay to T=0" mass fractions.

For mixed-element fits (e.g. Mn-56 + Fe shielding + …), the per-nuclide enrichments give the
mass fraction *within each chemical element*. For single-isotope elements (one Mn nuclide, one
Fe nuclide, …) that fraction is trivially 1.

The user typically configures the fit via:

- **ROIs**: energy ranges to fit, with a continuum form per ROI (linear, quadratic, stepped, …).
  See §3.7.
- **Rel-eff curve(s)**: see §1.2 below.
- **Energy-cal type**: `NoFit` / `LinearFit` / `NonLinearFit` (whether and how aggressively to
  adjust the spectrum's energy calibration during the fit).
- **FWHM form**: how peak widths scale with energy (`Gadras`, `Polynomial_2..6`, `Bernstein_2..6`).
- **FWHM estimation method**: where the initial FWHM comes from (DRF, peaks in spectrum, …).
- **Peak skew**: optional asymmetric-tail shape (`GaussExp`, `CrystalBall`, …).
- **Per-curve nuclide list**: nuclides + ages (and whether to fit each age) for that curve.

### 1.2 Relative-efficiency curve forms

A "relative-efficiency curve" `Eff(E)` describes the energy-dependent relative response of the
detector + intervening attenuation. The fit uses peak-area / yield / `Eff(E)` to recover the
relative activities; without an `Eff(E)` form there is no way to compare peak areas at different
energies. **One** rel-eff curve is the common case — it ties together every nuclide in the
fit and lets a single set of activity ratios explain every peak. **Multiple** curves are
allowed (in principle arbitrarily many; in practice no more than two), e.g. when the same
nuclide appears in two physically distinct populations that should be separated. A successful
solution always has at least one curve.

There are two families of rel-eff form, selected by `rel_eff_eqn_type`:

#### 1.2.1 Equation forms

A polynomial in some transform of energy. The form picks the transform; `rel_eff_eqn_order`
picks the polynomial order; `coefficients` are the fitted polynomial coefficients
(`a_0 + a_1·t + a_2·t² + …` where `t` is the transformed energy).

| `rel_eff_eqn_type` | Equation (`E` in keV) |
|---|---|
| `"LnX"`            | `Eff(E) = a_0 + a_1·ln(E) + a_2·ln(E)² + …`              |
| `"LnY"`            | `Eff(E) = exp(a_0 + a_1·E + a_2·E² + …)`                  |
| `"LnXLnY"`         | `Eff(E) = exp(a_0 + a_1·ln(E) + a_2·ln(E)² + …)`          |
| `"FRAM Empirical"` | FRAM-specific empirical form (5 parameters), see SAND-FRAM docs |

`equation_text` / `equation_html` / `js_rel_eff_eqn` are pre-built renderings of the fitted
curve in human-readable text, HTML (with Greek letters etc. encoded), and JavaScript-callable
form respectively. Templates rendering an interactive plot use `js_rel_eff_eqn` directly inside
a `<script>` tag (Inja does not auto-escape; this is intentional for these fields).

#### 1.2.2 Physical model (`"FRAM Physical"`)

A physically-motivated form built from:

- a **self-attenuator** (the source material itself, with finite areal density `μ·t`):
  `(1 − exp(−μ·t)) / (μ·t)`;
- one or more **external attenuators** (intervening shielding layers): `Π_i exp(−μ_i·t_i)`;
- the **detector intrinsic efficiency** `[Det. Eff.]` (from the DRF);
- an optional **Hoerl correction**: `E^b · c^(1/E)` with fitted `b`, `c`, that absorbs
  residual energy-dependent mismatch between the assumed physical model and the data.

For Physical-model curves, `coefficients` is a packed vector of the fitted parameters in the
order the solver uses internally; for templates, the human-readable representation is in
`equation_text` (e.g. `(1 - exp(-0*μ_U3O8))/(0*μ_U3O8) * [exp(-0*μ_Iron)] * [Det. Eff.] *
[E^0.0514 * 0.995^(1/E)]`). `js_rel_eff_eqn` is a JS function literal that evaluates the same
form. `rel_eff_eqn_order` is omitted for Physical-model curves.

The numeric attenuator names visible in `equation_text` (e.g. `μ_U3O8`, `μ_Iron`) come from
the user-supplied material names in the original fit setup. A template that wants to display
attenuator details should rely on `equation_text` / `equation_html` rather than try to extract
them from `coefficients`, as the coefficient packing is implementation-defined.

## 2. How a report is rendered

The C++ side ([`InterSpec/RelActAutoReport.h`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec/RelActAutoReport.h),
[`src/RelActAutoReport.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/RelActAutoReport.cpp))
is a three-step pipeline that each caller wires up:

```cpp
// 1. Build the JSON payload (full data; see §3 for shape).
const nlohmann::json data = RelActAutoReport::solution_to_json(solution);

// 2. Build an inja::Environment with the bundled templates and callbacks pre-registered.
//    `include_dir` is the search root for any {% include %} directives in user templates.
inja::Environment env = RelActAutoReport::get_default_inja_env(include_dir);

// 3. Render the chosen template against `data` using `env`.
RelActAutoReport::render_template(out_stream, env, data, "html", include_dir);
```

`get_default_inja_env(include_dir)` enables `set_trim_blocks(true)`, registers all the custom
callbacks (see §6), and pre-registers the two default templates as named includes:
- `default-rel-act-auto-html-results`
- `default-rel-act-auto-txt-results`

`render_template` reads the `tmplt` parameter to decide what to render:

```cpp
// Built-ins:
RelActAutoReport::render_template(out, env, data, "html");   // default HTML report
RelActAutoReport::render_template(out, env, data, "txt");    // built-in text
RelActAutoReport::render_template(out, env, data, "json");   // raw JSON dump (skips Inja)

// User template:
RelActAutoReport::render_template(out, env, data,
                                  "my_summary.tmplt.csv",    // filename or absolute path
                                  "/path/to/templates");     // include_dir (must match step 2)
```

Pulling JSON construction, env construction, and rendering apart lets callers (a) build the
JSON once and render multiple templates against it, and (b) extend `data` with extra keys
before rendering.

The `tmplt` parameter is the only knob that picks the template; accepted values:

| `tmplt` | What you get |
|---|---|
| `""` (the default) | the built-in HTML report (same as `"html"`) |
| `"html"` | `std_rel_eff_summary.tmplt.html` (per-file) |
| `"txt"` or `"text"` | `std_rel_eff_summary.tmplt.txt` (per-file) |
| `"json"` | raw JSON dump of the data (skips Inja entirely) |
| `"html-summary"` or `"summary"` | `std_multi_file_summary.tmplt.html` (multi-file batch summary; see §11) |
| any other string | treated as a path: looked up in `include_dir` first (if non-empty), then in `<InterSpec::writableDataDirectory()>/IsotopicsByNuclidesReportTmplts/`, then as an absolute path, then in the default template directory - if you are creating a new template file, it will be specified using this option |

The `include_dir` parameter is purely a search root for relative `tmplt` paths and for any
`{% include %}` directives in user templates. Defaults to `""` (no user-side includes).
`solution_to_json` always populates every field — there are no "what to include" toggles. See
§3.9 for which keys carry the heavy JS/CSS blobs and how to strip them when dumping JSON.

## 3. The JSON data model

The data passed to the template is documented by example. Every key listed below is *always*
present in a real JSON object used to render a template; §3.9 shows a complete payload
with the heaviest blobs elided for readability.

### 3.1 Top-level keys

| Key | Type | Notes |
|---|---|---|
| `status` | object | `success` (bool), `status_code` (int), `fail_reason` (string), `error_message` (string) |
| `spectrum_title` | string | The user-supplied analysis title |
| `have_multiple_rel_eff` | bool | True if more than one rel-eff curve was fit |
| `live_time_s` | number | Spectrum live time, seconds |
| `chi2`, `dof`, `chi2_per_dof` | number, int, number | Goodness of fit |
| `chi2_str`, `chi2_per_dof_str` | string | Pre-formatted with `%.6G` |
| `rel_eff_curves` | array | Per-curve equation info (see §3.2) |
| `rel_eff_chart_json` | string | JSON literal interpolated into `<script>` to drive the bundled `RelEffPlot` widget; HTML reports only |
| `relative_activities` | array | Per-curve activity rows; nested `pu` and `ratios` (see §3.3) |
| `pu_corrections` | array | Same Pu rows as `relative_activities[i].pu`; convenience for templates that iterate Pu data without nesting through `relative_activities` |
| `ratios` | array | Same ratio rows as `relative_activities[i].ratios`; convenience |
| `peaks` | array | Per-peak rows (see §3.4) |
| `rois` | array | ROI ranges (see §3.7) |
| `energy_calibration` | object | Adjustments + deviation pairs (see §3.5) |
| `options` | object | Top-level fit options (`fwhm_form`, `energy_cal_type`, `skew_type`, etc.) |
| `curve_options` | array | Per-curve fit options (`rel_eff_eqn_type`, `rel_eff_eqn_order`, `nucs_of_el_same_age`) |
| `foreground`, `background` | object | Spectrum metadata: `filename`, `title`, `sample_numbers` (array of int), `display_filename` (pre-formatted "Filename, Samples: {…}"), `live_time`, `real_time`, `start_time`. See note below. |
| `detector` | object | DRF metadata (name, description, diameter_mm, has_fwhm, has_intrinsic_eff) |
| `timing` | object | Solver timings + duration string |
| `warnings`, `warnings_html` | array of string | HTML-sanitized variant for `warnings_html` |
| `timestamps` | object | `local`, `utc` |
| `compile_timestamp` | string | `__TIMESTAMP__` of the executable |
| `assets` | object | `D3_JS`, `SpectrumChart_JS`, `SpectrumChart_CSS`, `RelEffPlot_JS`, `RelEffPlot_CSS` — multi-hundred-KB JS/CSS strings, only useful for HTML reports |
| `spectrum_chart` | object | `have_spectrum` (bool); when true also `div_id`, `set_js`, `lower_energy`, `upper_energy`. HTML reports only. |

### 3.2 `rel_eff_curves[i]`

A successful solution always has **at least one** entry. The single-curve case dominates by
far; two-curve fits are rare. Templates that don't care about multi-curve presentation can
just iterate the array and ignore `have_multiple_rel_eff`. Templates that *do* care typically
gate per-curve headings on `have_multiple_rel_eff`.

Each entry has these fields (see §1.2 for what the equation forms mean):

```json
{
  "index": 0,
  "name": "Curve 0",
  "rel_eff_eqn_type": "LnX",            // "LnX" | "LnY" | "LnXLnY" | "FRAM Empirical" | "FRAM Physical"
  "rel_eff_eqn_order": 3,                // omitted for FRAM Physical
  "coefficients": [692.93, -284.58, 38.96, -1.78],
  "equation_text": "692.93 - 284.58*ln(x)^1 + 38.96*ln(x)^2 - 1.7756*ln(x)^3",
  "equation_html": "692.93 - 284.58*ln(x)^1 + 38.96*ln(x)^2 - 1.7756*ln(x)^3",
  "js_rel_eff_eqn": "function(x){ return 692.928444 + -284.582067*Math.log(x) + ... ; }"
}
```

| Field | Meaning |
|---|---|
| `index` | Zero-based index into the per-curve arrays (`relative_activities`, `peaks[i].rel_eff_index`, etc.) |
| `name` | User-supplied label for this curve (defaults to `"Curve N"`) |
| `rel_eff_eqn_type` | See §1.2 — selects the equation family or Physical model |
| `rel_eff_eqn_order` | Polynomial order; omitted entirely for `FRAM Physical` |
| `coefficients` | Fitted polynomial coefficients (Equation forms) or packed Physical-model parameters (treat as opaque) |
| `equation_text` | Plain-text rendering of the fitted curve, suitable for monospace text reports |
| `equation_html` | HTML-escaped rendering, with Greek letters / sub/superscript markup as appropriate |
| `js_rel_eff_eqn` | JavaScript function literal — interpolate **raw** into a `<script>` tag (Inja does not auto-escape; this is intentional) |

For Physical-model curves, `equation_text` / `equation_html` are the only human-readable view
of the fitted self/external attenuators and Hoerl coefficients (see §1.2.2). The numeric
parameters are packed into `coefficients` in an internal order; do not try to parse them in a
template.

### 3.3 `relative_activities[i]` (per curve)

One entry per rel-eff curve. Each entry holds the per-nuclide fit results for that curve, plus
a nested `pu` block (Pu fits only) and a nested `ratios` block (when the curve has ≥2 nuclides).

```json
{
  "curve_index": 0,
  "nuclides": [
    {
      "name": "Mn56",
      "rel_activity": 14.37,
      "rel_activity_uncertainty": 0.55,
      "rel_activity_with_uncert_str": "14.3692 ± 0.552453",
      "age": 0.0, "age_uncertainty": 0.0, "age_was_fit": false,
      "age_str": "0",
      "is_pu242_by_corr": false,
      "rel_mass": 1.79e-17,
      "total_mass_fraction": 1.0,        // 0..1, multiply by 100 for %
      "enrichment": 1.0,                  // 0..1
      "has_enrichment_uncert": true,
      "enrichment_uncert": 0.0,
      "enrichment_minus_2sigma": 1.0,
      "enrichment_plus_2sigma": 1.0,
      "detector_counts": 36847.4
    }
  ],
  "pu":     { /* same shape as pu_corrections[i] - see §3.8 - present only for Pu fits */ },
  "ratios": { /* same shape as ratios[i] - see below - present iff there are >=2 nuclides */ }
}
```

#### Nuclide entry fields

| Field | Meaning |
|---|---|
| `name` | Nuclide symbol (`"U235"`, `"Mn56"`) — or element / reaction name for non-nuclide sources |
| `rel_activity` | The fitted activity of this nuclide in arbitrary units; only ratios between nuclides are physically meaningful |
| `rel_activity_uncertainty` | 1σ uncertainty on `rel_activity` (in the same arbitrary units) |
| `rel_activity_with_uncert_str` | Pre-formatted "value ± uncert" string |
| `age` | Fitted (or fixed) source age in seconds. Used by SandiaDecay to age daughter contributions |
| `age_uncertainty` | 1σ uncertainty on `age` (seconds), only meaningful when `age_was_fit` is true |
| `age_was_fit` | True if the solver fit this nuclide's age, false if held fixed |
| `age_str` | Pre-formatted age string in best-fitting time units (e.g. `"20 years"`, `"3.2 days"`, `"0"`) |
| `is_pu242_by_corr` | True only for the Pu242 row when Pu242 was estimated by correlation rather than directly fit |
| `rel_mass` | `rel_activity / activity_per_gram(nuclide)` — relative mass in arbitrary units |
| `total_mass_fraction` | This nuclide's mass divided by the sum of all nuclide masses in the curve. 0..1; multiply by 100 for %. (Often less informative than `enrichment` for mixed-element fits.) |
| `enrichment` | This nuclide's mass divided by the total mass of *its element* across the fit (this is the "enrichment" or per-element mass fraction). 0..1; multiply by 100 for %. For a single-isotope element this is trivially 1 |
| `has_enrichment_uncert` | True when covariance was available and a 1σ enrichment uncertainty could be computed |
| `enrichment_uncert` | 1σ enrichment uncertainty (0..1 scale) — only meaningful when `has_enrichment_uncert` |
| `enrichment_minus_2sigma`, `enrichment_plus_2sigma` | 2σ enrichment band (0..1); only meaningful when `has_enrichment_uncert` |
| `detector_counts` | Total fitted detector counts attributed to this nuclide across all its peaks in this curve |

`enrichment` is what U/Pu reports usually want to display ("99.3% U-238", "0.66% U-235"), and
is also the right field to drive the §1.1 "per-element mass-fraction breakdown" idea — when
the value is essentially 1.0 (e.g. > 0.999), the nuclide is the only one of its element in the
fit and the breakdown is trivially uninteresting.

#### `ratios[i]` (and the nested `relative_activities[i].ratios`)

Present iff a curve has ≥ 2 nuclides. Lists *every unique unordered pair* of nuclides in that
curve. Each pair carries both directions of the activity ratio and mass ratio, plus their
uncertainties — *both directions are pre-computed server-side*, so templates never have to
divide.

```json
{
  "curve_index": 0,
  "pairs": [
    {
      "numerator": "U238", "denominator": "U235",
      "activity_ratio": 23.40, "activity_ratio_str": "23.4021",
      "activity_ratio_inv": 0.0427, "activity_ratio_inv_str": "0.0427311",
      "mass_ratio": 150.42, "mass_ratio_str": "150.423",
      "mass_ratio_inv": 0.00665, "mass_ratio_inv_str": "0.00664791",
      "activity_ratio_uncertainty": 0.001, "activity_ratio_uncertainty_inv": 0.001,
      "uncert_pct_str": "0.123%", "uncert_inv_pct_str": "0.123%"
    }
  ]
}
```

If covariance was unavailable, `uncert_pct_str` and `uncert_inv_pct_str` are `"--"`.

### 3.4 `peaks[i]`

Sorted by energy, ascending. Each entry is one fit Gaussian (plus its continuum) in some ROI.
Two flavours appear, distinguished by `is_floating`:

- **Nuclide-assigned peaks** (`is_floating: false`) — the fit constrained the peak's amplitude
  through the relative-activity / rel-eff machinery. `nuclide` names the source (a nuclide
  symbol, an element symbol for x-rays, or a reaction name); `yield`, `cps_over_yield`,
  `fit_rel_eff`, and `fit_rel_eff_uncert_pct` are meaningful.
- **Free-floating peaks** (`is_floating: true`) — peaks the user explicitly told the solver to
  fit *without* tying their amplitude to any nuclide's relative activity. They appear in the
  spectrum (they were observed) and need to be modelled so they don't pollute neighbouring
  ROIs, but their amplitude is determined directly from the data alone. The user typically
  adds these for genuinely unidentified peaks, single-escape / double-escape peaks they don't
  want to attribute, or contaminants. For floating rows: `nuclide` is `""`, `yield` and
  `cps_over_yield` are 0, `fit_rel_eff` and `fit_rel_eff_uncert_pct` are -1, only the energy
  / amplitude / continuum columns are meaningful.

Templates rendering a peak table typically branch on `is_floating` to decide which columns to
fill (see e.g. `std_rel_eff_summary.tmplt.html`).

```json
{
  "energy": 846.75,
  "amplitude": 26178.91,
  "amplitude_uncertainty": 164.80,
  "amplitude_uncertainty_percent": 0.63,    // already multiplied by 100
  "yield": 0.989,                            // 0 for floating peaks
  "cps_over_yield": 14.71,                   // 0 for floating peaks
  "fit_rel_eff": 1.0245,                     // -1 for floating peaks or if covariance unavailable
  "fit_rel_eff_uncert_pct": 1.68,            // -1 if uncertainty unavailable
  "rel_eff_index": 0,
  "is_floating": false,
  "nuclide": "Mn56",                         // or element symbol or reaction name; "" for floating
  "source_type": "nuclide",                  // "nuclide" | "element" | "reaction" | "floating"
  "continuum": {
    "type": "pct-linear",                    // "Linear" | "Quadratic" | "Cubic" | "FlatStep" | "LinearStep" | "BiLinearStep" | "External" | "NoOffset"
    "lower_energy": 836.45,
    "upper_energy": 855.05
  }
}
```

| Field | Meaning |
|---|---|
| `energy` | Fitted peak mean (keV) |
| `amplitude` | Fitted peak area (counts) |
| `amplitude_uncertainty` | 1σ uncertainty on the area (counts) |
| `amplitude_uncertainty_percent` | `100 · amplitude_uncertainty / amplitude` — already a percentage |
| `yield` | Branching ratio × age-corrected intensity for this gamma; 0 for floating peaks |
| `cps_over_yield` | `amplitude / (yield · live_time)` — useful for visual rel-eff inspection; 0 for floating peaks |
| `fit_rel_eff` | Value of the fitted rel-eff curve at this peak's energy; -1 for floating or if covariance was unavailable |
| `fit_rel_eff_uncert_pct` | Percent 1σ uncertainty on `fit_rel_eff`; -1 if unavailable |
| `rel_eff_index` | Which `rel_eff_curves[i]` curve this peak belongs to |
| `is_floating` | See above |
| `nuclide`, `source_type` | Identification of the source; `source_type` is `"floating"` for free-floaters |
| `continuum.type` | Continuum form fit beneath the peak in this ROI (see §3.7 for the full list of values and what they mean) |
| `continuum.lower_energy`, `continuum.upper_energy` | The ROI bounds *as the solver actually used them* (after any FWHM-based expansion) |

### 3.5 `energy_calibration`

```json
{
  "was_fit": true,
  "energy_cal_type": "NonLinearFit",
  "adjustments": [
    {
      "parameter_index": 1,                       // 0=offset, 1=gain, 2=cubic
      "physical_value": 2.228e-05,
      "physical_value_str": "2.22803E-05",
      "raw_value": 1.0036504066868457,
      "type": "gain",                              // "offset" | "gain" | "cubic"
      "type_long": "gain adjustment",
      "units": "keV/chnl",                         // short
      "units_long": "keV/channel",                 // long (HTML-friendly)
      "was_fit": true
    }
  ],
  "deviation_pair_offsets": [
    { "anchor": 60.0, "offset": 0.0, "display": "{60.0,0.00}" },
    { "anchor": 356.0, "offset": -0.59, "display": "{356.0,-0.59}" }
  ]
}
```

### 3.6 `foreground` and `background`

Spectrum-file metadata identifying *where the data came from*. `foreground` is always present
when a spectrum was passed to the fit; `background` is present only if the user supplied a
background spectrum.

```json
"foreground": {
  "filename":          "field_meas.n42",
  "title":             "Foreground",
  "sample_numbers":    [1, 3, 5],
  "display_filename":  "field_meas.n42, Samples: {1,3,5}",
  "live_time":         1800.0,
  "real_time":         1824.65,
  "start_time":        "20131226T211500"
}
```

(The example above shows the non-trivial case where the foreground was summed from a subset of
samples — `sample_numbers` is populated and `display_filename` includes the `Samples: {…}`
annotation. For a single-sample file or a file whose fg+bg cover every sample, `sample_numbers`
would be `[]` and `display_filename` would just be the filename — see the rule below.)

| Field | Meaning |
|---|---|
| `filename` | Basename (no path) of the file the spectrum was loaded from. Empty when the caller didn't supply one |
| `title` | The user-supplied title on the foreground/background `Measurement` (or the file's title if the user didn't set one) |
| `sample_numbers` | The sample numbers within `filename` that were summed to produce the foreground/background. Empty when the sample structure is trivial — see below |
| `display_filename` | Pre-formatted convenience string: just `filename` when `sample_numbers` is empty, otherwise `"filename, Samples: {1,3,5}"`. Templates that just want a label use this |
| `live_time`, `real_time` | Detector live/real time, seconds |
| `start_time` | ISO-8601 acquisition start time (foreground only) |

#### When are `sample_numbers` empty?

`sample_numbers` is left empty (and the `display_filename` falls back to just the filename)
when the sample-number annotation would be redundant:

- the foreground (or background) uses *every* sample in its file (most commonly: a single-sample
  file); or
- the foreground and background read from the same file and *together* cover every sample in
  it (the typical "two-sample file: sample 1 = fg, sample 2 = bg" case).

In every other case `sample_numbers` lists the samples actually used. Templates that want to
unconditionally show the samples can iterate `sample_numbers` themselves; templates that want
the default "only show when non-trivial" behaviour can just use `display_filename`.

### 3.7 `rois[i]`

A "Region of Interest" is a contiguous energy range the user marked for fitting. Each ROI
contains one or more peaks; the solver fits the Gaussian peak parameters on top of a continuum
of the chosen form. ROIs are user-supplied (typically by selecting peaks on the spectrum); the
fit does not invent ROIs.

```json
{
  "lower_energy": 836.4,
  "upper_energy": 855.0,
  "continuum_type": "pct-linear",
  "range_limits_type": "Fixed"
}
```

| Field | Meaning |
|---|---|
| `lower_energy`, `upper_energy` | ROI energy bounds (keV) as configured by the user, *before* any FWHM-based expansion the solver may have applied |
| `continuum_type` | Continuum form fit beneath the peaks in this ROI. See the table below |
| `range_limits_type` | How the solver treats the ROI bounds during the fit — see below |

#### `continuum_type` values

A continuum is the slowly-varying background under the Gaussian peaks in the ROI. The same
strings appear in `peaks[i].continuum.type`.

| Value | What it means |
|---|---|
| `"NoOffset"`     | No continuum (peaks sit on zero) |
| `"Constant"`     | Flat continuum (`a_0`) |
| `"Linear"`       | Linear in energy (`a_0 + a_1·E`) |
| `"Quadratic"`    | Quadratic in energy |
| `"Cubic"`        | Cubic in energy |
| `"FlatStep"`     | Flat with a step at peak mean — for high-statistics peaks where the continuum drops across each peak |
| `"LinearStep"`   | Linear with a step |
| `"BiLinearStep"` | Two linear segments meeting at the peak mean |
| `"External"`     | Continuum determined by an external fit (rare) |

The string `"pct-linear"` and similar dashed/lower-case forms appear when the spectrum file's
internal i18n key was not translated (the solver was run outside a Wt session); they are
synonymous with the canonical capitalised form for that continuum type.

#### `range_limits_type` values

| Value | What it means |
|---|---|
| `"Fixed"`            | Bounds are strictly enforced — the solver will not move them |
| `"CanExpandForFwhm"` | Bounds may expand outward by a fraction of the FWHM at the edge if a peak's tail leaks past them |
| `"CanBeBrokenUp"`    | The ROI may be subdivided into multiple smaller ROIs based on where peaks actually fall |

### 3.8 `pu_corrections[i]` (only when Pu nuclides present)

```json
{
  "curve_index": 0,
  "rows": [
    { "nuclide": "Pu239", "mass_number": 239, "is_by_corr": false,
      "uncorrected_mass_frac": 0.93, "mass_frac": 0.92, "mass_frac_uncert": 0.001 },
    { "nuclide": "Pu242", "mass_number": 242, "is_by_corr": true,
      "uncorrected_mass_frac": -1, "mass_frac": 0.012, "mass_frac_uncert": 0.0006 }
  ],
  "mass_fractions": { "pu238":..., "pu239":..., "pu240":..., "pu241":..., "pu242":..., "pu242_uncert":..., "is_within_range": true },
  "uncorrected_mass_fractions": { "pu238":..., "pu239":..., "pu240":..., "pu241":..., "pu_other":..., "age":"...", "age_seconds":... },
  "back_decayed": [ { "nuclide": "Pu239", "mass_frac": 0.94 }, ... ]   // present only when ages are uniform and >1 Pu nuclide
}
```

`uncorrected_mass_frac == -1` and `mass_frac_uncert < 0` are both blank-cell sentinels.

### 3.9 Full Mn56 example

> **Note — heavy blobs elided.** `solution_to_json` *always* populates the JS/CSS asset blobs
> and the spectrum-chart init script; the example below has had these stripped for readability.
> A real `tmplt = "json"` dump additionally contains the following keys, each holding a
> multi-hundred-KB string of literal JS/CSS:
>
> - `assets.D3_JS`
> - `assets.SpectrumChart_JS`
> - `assets.SpectrumChart_CSS`
> - `assets.RelEffPlot_JS`
> - `assets.RelEffPlot_CSS`
> - `spectrum_chart.set_js`
>
> The remaining `spectrum_chart.*` keys (`have_spectrum`, `div_id`, `lower_energy`,
> `upper_energy`) are small and worth keeping. `rel_eff_chart_json` is a JSON-encoded string
> of the rel-eff chart payload — typically a few KB, kept by default.
>
> To produce a clean dump like the one below, pipe through `jq`:
>
> ```bash
> jq 'del(.assets, .spectrum_chart.set_js)' mn56_inja.json > mn56_inja.clean.json
> ```

```json
{
  "assets": {},
  "chi2": 167.355, "chi2_per_dof": 0.764, "chi2_str": "167.355", "chi2_per_dof_str": "0.764176",
  "compile_timestamp": "Sat Apr 25 20:45:31 2026",
  "curve_options": [
    { "name": "Curve 0", "nucs_of_el_same_age": false, "rel_eff_eqn_order": 3, "rel_eff_eqn_type": "LnX" }
  ],
  "detector": {
    "description": "65 x 50-mm, SimpleMass - from Relative Eff. File",
    "diameter_mm": 65.0, "has_fwhm": true, "has_intrinsic_eff": true,
    "name": "ORTEC Detective-X_LANL_100cm (59%)"
  },
  "dof": 219,
  "energy_calibration": {
    "adjustments": [
      { "parameter_index": 1, "physical_value": 2.228e-05, "physical_value_str": "2.22803E-05",
        "raw_value": 1.0036504, "type": "gain", "type_long": "gain adjustment",
        "units": "keV/chnl", "units_long": "keV/channel", "was_fit": true }
    ],
    "deviation_pair_offsets": [
      { "anchor": 60.0,  "display": "{60.0,0.00}",  "offset": 0.0 },
      { "anchor": 356.0, "display": "{356.0,-0.59}", "offset": -0.59 },
      { "anchor": 583.0, "display": "{583.0,-0.75}", "offset": -0.75 }
    ],
    "energy_cal_type": "NonLinearFit",
    "was_fit": true
  },
  "foreground": {
    "filename": "Mn56 Shielded", "live_time": 1800.0, "real_time": 1824.65,
    "start_time": "20131226T211500", "title": "Mn56 Shielded"
  },
  "have_multiple_rel_eff": false,
  "live_time_s": 1800.0,
  "options": {
    "additional_br_uncert": 0.05, "energy_cal_type": "NonLinearFit",
    "fwhm_estimation_method": "StartFromDetEffOrPeaksInSpectrum",
    "fwhm_form": "Polynomial_4", "skew_type": 2, "skew_type_str": "GaussExp",
    "spectrum_title": "Foreground"
  },
  "peaks": [
    {
      "amplitude": 26178.9, "amplitude_uncertainty": 164.8, "amplitude_uncertainty_percent": 0.63,
      "continuum": { "lower_energy": 836.45, "type": "pct-linear", "upper_energy": 855.05 },
      "cps_over_yield": 14.71, "energy": 846.75, "fit_rel_eff": 1.0245, "fit_rel_eff_uncert_pct": 1.68,
      "is_floating": false, "nuclide": "Mn56", "rel_eff_index": 0, "source_type": "nuclide",
      "yield": 0.989
    }
  ],
  "pu_corrections": [],
  "ratios": [ { "curve_index": 0, "pairs": [] } ],
  "rel_eff_chart_json": "[]",
  "rel_eff_curves": [
    {
      "coefficients": [692.93, -284.58, 38.96, -1.776],
      "equation_html": "692.93 - 284.58*ln(x)^1 + 38.96*ln(x)^2 - 1.7756*ln(x)^3",
      "equation_text": "692.93 - 284.58*ln(x)^1 + 38.96*ln(x)^2 - 1.7756*ln(x)^3",
      "index": 0,
      "js_rel_eff_eqn": "function(x){ return 692.93 + -284.58*Math.log(x) + ... ; }",
      "name": "Curve 0", "rel_eff_eqn_order": 3, "rel_eff_eqn_type": "LnX"
    }
  ],
  "relative_activities": [
    {
      "curve_index": 0,
      "nuclides": [
        {
          "age": 0.0, "age_str": "0", "age_uncertainty": 0.0, "age_was_fit": false,
          "detector_counts": 36847.4, "enrichment": 1.0,
          "enrichment_minus_2sigma": 1.0, "enrichment_plus_2sigma": 1.0, "enrichment_uncert": 0.0,
          "has_enrichment_uncert": true, "is_pu242_by_corr": false, "name": "Mn56",
          "rel_activity": 14.37, "rel_activity_uncertainty": 0.55,
          "rel_activity_with_uncert_str": "14.3692 ± 0.552453",
          "rel_mass": 1.79e-17, "total_mass_fraction": 1.0
        }
      ],
      "ratios": { "curve_index": 0, "pairs": [] }
    }
  ],
  "rois": [
    { "continuum_type": "pct-linear", "lower_energy": 836.4, "range_limits_type": "Fixed",
      "upper_energy": 855.0 }
  ],
  "spectrum_chart": { "have_spectrum": false },
  "spectrum_title": "Foreground",
  "status": { "error_message": "", "fail_reason": "Success", "status_code": 0, "success": true },
  "timestamps": { "local": "26-Apr-2026 08:53:25 PM", "utc": "27-Apr-2026 03:53:25" },
  "timing": {
    "cov_calls": 2, "duration_str": "87.76 ms", "function_eval_solution": 197,
    "function_eval_total": 199, "microseconds_eval": 87759, "microseconds_in_eval": 21826,
    "solve_calls": 197
  },
  "warnings": [], "warnings_html": []
}
```

Generate a fresh, identical-shape JSON for any solution by passing `tmplt = "json"` to
`render_template` (or just call `data.dump(2)` on the JSON returned by `solution_to_json`),
then strip the heavy blobs with the `jq` command above.

## 4. Inja syntax — what works in this build

InterSpec ships Inja 3.x (in `external_libs/SpecUtils/3rdparty/inja/inja.hpp`). Reliable features:

```jinja
{# This is a comment #}

{{ variable }}                              {# interpolate #}
{{ object.field }}                          {# dotted access #}
{{ array.0 }}                               {# numeric-literal indexing into an array #}
{{ at(array, n) }}                          {# variable indexing - use this for `array[expr]` #}

{% if status.success %} ok {% endif %}
{% if foo %} ... {% else if bar %} ... {% else %} ... {% endif %}

{% for item in collection %}
    {{ item.name }}
    {% if loop.is_first %}-- first --{% endif %}
    {% if loop.is_last  %}-- last  --{% endif %}
    index = {{ loop.index }}                {# 0-based #}
    {{ at(parent_array, loop.index).field }}{# variable index via at() #}
{% endfor %}

{% include "default-rel-act-auto-html-results" %}        {# named, pre-registered #}

{# Boolean-y predicates #}
{% if length(warnings) > 0 %}
{% if existsIn(curve_acts, "pu") %}
{% if status.status_code == 0 %}

{# Numeric comparisons #}
{% if mass_frac_uncert >= 0 %} ... {% else %} ... {% endif %}
```

`set_trim_blocks(true)` is set, so a `{% ... %}` block on its own line removes the newline that
follows it. To preserve a newline, use a non-block tag (e.g. wrap in `{{ "" }}`).

## 5. What does NOT work (or is fragile)

These pitfalls have bitten us. Avoid them or work around them in `solution_to_json` rather than
in templates:

1. **Bracket array indexing (`arr[expr]`) doesn't work — use `at(arr, expr)` instead.**
   `Files[loop.index]` triggers a JSON-parse error in Inja's expression evaluator, but the
   built-in `at(array, index)` callback is fully supported (it's the documented Inja
   variable-indexing accessor — see [pantor/inja](https://pantor.github.io/inja/)). For
   numeric-literal indices, the dot form `arr.0` is shorter; for variables (`loop.index`,
   the result of another callback, etc.) use `at`:
   ```jinja
   {{ at(Files, loop.index).Filename }}            {# fine #}
   {{ at(curve_acts.nuclides, 0).name }}           {# also fine; equivalent to .nuclides.0.name #}
   {{ Files[loop.index].Filename }}                {# DOES NOT WORK -- JSON parse error #}
   ```
   Cross-referencing through nested C++ JSON (e.g. `relative_activities[i].pu` already nests
   the per-curve Pu data inside the per-curve relative-activities object) still avoids needing
   the index altogether and is the cleaner choice when feasible.
2. **Arithmetic in expressions can be silently dropped or wrong.** `nuc.total_mass_fraction * 100`
   has rendered as `1` instead of `100` in some contexts. Prefer the `pct(x)` callback (§6) or
   pre-compute a `_str` field server-side.
3. **No HTML auto-escaping.** Strings interpolated with `{{ }}` are emitted verbatim. Use
   `{{ safe_html(s) }}` for any string that might contain `<`, `>`, `&`, or quotes (e.g.
   warnings or user-supplied titles). Embedded JS/CSS blobs like `assets.SpectrumChart_JS` are
   *meant* to be raw — interpolate them inside `<script>` / `<style>` tags as-is.
4. **Compound `{% if a and b %}`** can have ordering oddities; prefer nested `{% if %}` blocks.
5. **`{% include %}` requires a registered name or, with a custom template directory, a
   filename relative to that directory.** Web builds (`BUILD_FOR_WEB_DEPLOYMENT`) disable
   filesystem includes via `set_search_included_templates_in_files(false)`.
6. **Template-side division.** Reciprocals like `1.0 / pair.activity_ratio` work in trivial cases
   but are risky for zero/negative divisors. Both directions of every ratio are pre-computed
   server-side (`activity_ratio_inv`, `mass_ratio_inv`, `*_inv_str`).
7. **`set_trim_blocks(true)` is aggressive — eats whitespace, not just newlines.** This Inja
   version's trim-blocks consumes *every* space and tab immediately after a `{% ... %}` close,
   then optionally one newline. So `{% endif %} | next` becomes `{% endif %}|next` (the leading
   space before `|` is gone), and `{% endif %}\n` at end of line eats the newline (the next line
   joins the current one). Workarounds:
   - **End-of-line `{% endif %}` /  `{% endfor %}`** — append an empty interpolation `{{ "" }}`
     after the close: `{% endif %}{{ "" }}\n`. Inja's trim only fires on `{% %}` closes, not on
     `{{ }}` closes, so the trailing newline is preserved.
   - **Mid-line `{% endif %}` followed by literal whitespace + content** — move the trailing
     whitespace *inside* the if/else branches:
     `{% if cond %}A    {% else %}B    {% endif %}NEXT` instead of
     `{% if cond %}A{% else %}B{% endif %}    NEXT` (the four spaces between `{% endif %}` and
     `NEXT` would be eaten by trim).
   - **Markdown headings** beware: `## ` and `### ` at column 0 are *also* the Inja line-statement
     prefix, so `## Section` parses as a malformed Inja statement. Use `<h2>...</h2>` HTML inside
     Markdown templates instead (Markdown allows raw HTML).

## 6. Custom Inja callbacks (registered by `get_default_inja_env`)

These are the callable helpers available inside `{{ ... }}` and `{% ... %}` expressions.
**Use them instead of in-template arithmetic** — see §5 for why expression arithmetic in
this Inja version is fragile. `existsIn` and `length` are inja built-ins; the rest are
InterSpec-specific.

| Callback | Args | Returns | Example |
|---|---|---|---|
| `at(array, index)` | `array, int` | element | `at(Files, loop.index).Filename` (variable indexing — see §5.1) |
| `printFixed(value, decimals)` | `number, int` | string | `printFixed(p.energy, 2)` -> `"846.75"` |
| `printCompact(value, sigfigs)` | `number, int` | string | `printCompact(amp, 6)` -> `"26178.9"` |
| `pct(value)` or `pct(value, sigfigs)` | `number [, int]` | string | `pct(nuc.enrichment, 4)` -> `"99.33"` |
| `safe_html(str)` | `string` | string | `safe_html(warnings.0)` |
| `scientific(value, prec)` | `number, int` | string | `scientific(rel_act, 6)` -> `"1.43692e+01"` |
| `format(fmt, value)` | `string, number` | string | `format("%.4f", chi2)` |
| `default(val, fallback)` | `value, value` | string | `default(opts.title, "Untitled")` |
| `existsIn(obj, "key")` | `object, string` | bool | `existsIn(nuc, "enrichment")` |
| `length(arr)` | `array or string` | int | `length(warnings)` |

### How to use them

- **`at(arr, idx)`** — array indexing with a *variable* index (e.g. `at(Files, loop.index)`,
  `at(arr, idx + 1)`). Bracket form `arr[expr]` is not supported in expression position
  (see §5.1). For numeric-literal indices the `arr.0` form is shorter; use `at()` whenever the
  index is anything other than a constant. Returns the element by reference, so chaining is
  fine: `{{ at(Files, loop.index).relative_activities.0.nuclides.0.name }}`.

- **`printFixed(x, n)`** — fixed-point with `n` decimals. Use for angles, energies (keV), and
  any quantity where the column needs a stable number of decimals. `printFixed(846.7540283, 2)`
  → `"846.75"`. For null / NaN / infinite inputs returns `"--"`.

- **`printCompact(x, n)`** — `n` significant figures, choosing between fixed and scientific
  notation automatically. The general-purpose numeric formatter. `printCompact(26178.91, 6)`
  → `"26178.9"`; `printCompact(6.67e-10, 4)` → `"6.67E-10"`. For null / NaN / infinite returns
  `"--"`. **Almost always the right choice when displaying a fitted value.**

- **`pct(x)` / `pct(x, n)`** — multiplies by 100 and formats with `printCompact` at `n`
  significant figures (default 4). Used because expression-side multiplication (`x * 100`)
  is unreliable in this Inja version. `pct(0.99333, 4)` → `"99.33"`. Append `%` yourself in
  the surrounding template (`{{ pct(nuc.enrichment, 4) }}%`).

- **`safe_html(s)`** — escapes `&`, `<`, `>`, `'`, `"` for safe insertion into HTML body
  text. **Required** for any string from `warnings`, user-supplied `spectrum_title`, or
  similar — Inja does not auto-escape. *Do not* apply to `assets.SpectrumChart_JS`,
  `assets.RelEffPlot_JS`, `js_rel_eff_eqn`, `spectrum_chart.set_js`, or `rel_eff_chart_json`
  — those are meant to interpolate raw inside `<script>` / `<style>` tags. Note that
  `warnings_html` is the same data as `warnings`, *already sanitized* by the C++ side; in
  modern templates prefer iterating `warnings_html` directly over calling `safe_html` on
  `warnings`.

- **`scientific(x, n)`** — forces scientific notation at `n` significant figures. Rarely
  needed in practice; `printCompact` chooses scientific automatically when values get small or
  large.

- **`format(fmt, x)`** — printf-style. Useful for unusual specifiers (`format("%5.2e", chi2)`)
  but **the format string is passed straight to `snprintf`**: only ever pass a literal string,
  never user-controlled data, or you create a format-string injection vulnerability.

- **`default(v, fallback)`** — if `v` is null / empty, returns `fallback`. Convenient for
  optional metadata fields. Most numeric callbacks already return `"--"` on null, so this is
  mostly useful for strings.

- **`existsIn(obj, "key")`** — true iff `obj` has the named key. Use to gate optional fields:
  `{% if existsIn(nuc, "age_str") %}{{ nuc.age_str }}{% endif %}`.

- **`length(x)`** — array length or string length. Use for "is the array empty" tests
  (`{% if length(warnings) > 0 %}`); inja's truthiness rules around empty arrays are
  unreliable, so this explicit form is recommended.

For null / NaN / infinite inputs, all numeric formatters (`printFixed`, `printCompact`, `pct`)
return the string `"--"` rather than throwing — this lets templates assume the call always
produces a printable string.

## 7. Common patterns

### 7.1 Iterating per-curve activities + Pu correction + ratios in one loop

```jinja
{% for curve_acts in relative_activities %}
  <h3>Curve {{ curve_acts.curve_index }}</h3>

  {% if existsIn(curve_acts, "pu") %}
    {# render Pu mass-fraction table from curve_acts.pu.rows #}
  {% endif %}

  <table>
    <tr><th>Nuclide</th><th>Rel. Act.</th><th>Mass %</th><th>Enrichment</th></tr>
    {% for nuc in curve_acts.nuclides %}
      <tr>
        <td>{{ nuc.name }}{% if nuc.is_pu242_by_corr %} (by corr){% endif %}</td>
        <td>{{ printCompact(nuc.rel_activity, 6) }} &plusmn; {{ printCompact(nuc.rel_activity_uncertainty, 6) }}</td>
        <td>{{ pct(nuc.total_mass_fraction, 4) }}%</td>
        <td>{{ pct(nuc.enrichment, 4) }}%</td>
      </tr>
    {% endfor %}
  </table>

  {% if existsIn(curve_acts, "ratios") %}
    {% if length(curve_acts.ratios.pairs) > 0 %}
      {# render ratio table #}
    {% endif %}
  {% endif %}
{% endfor %}
```

### 7.2 Embedding the rel-eff plot

> **HTML reports only.** This pattern produces an interactive D3 chart in a browser; it has no
> meaning in plain-text, CSV, or Markdown reports — those should consume `rel_eff_curves[i]`
> directly (e.g. show `equation_text`).

`rel_eff_chart_json`, `assets.RelEffPlot_JS`, and `assets.RelEffPlot_CSS` are always present in
the JSON; the template just needs to drop them into the right tags. The values are raw
JS/CSS/JSON literals — interpolate them with `{{ … }}` inside `<script>` / `<style>` and do
**not** pass them through `safe_html`.

```jinja
<style>{{ assets.RelEffPlot_CSS }}</style>
<script>{{ assets.RelEffPlot_JS }}</script>

<div id="releffchart" style="height: 30vw;"></div>
<script>
  const rel_eff_data = {{ rel_eff_chart_json }};
  const chart = new RelEffPlot("releffchart");
  chart.setRelEffData(rel_eff_data);
</script>
```

### 7.3 Embedding the spectrum chart

> **HTML reports only.** Same caveat as §7.2 — non-HTML reports should ignore the
> `spectrum_chart` and `assets.SpectrumChart_*` keys entirely.

The C++ side only populates the chart payload when a spectrum was provided to the fit. Always
check `spectrum_chart.have_spectrum` before emitting the chart markup; otherwise omit it (or
substitute a placeholder spacer).

The template writes its own `<div>` (with whatever layout / styling / class it likes); the
only constraint is the element's `id` must match `spectrum_chart.div_id`, because the bundled
`spectrum_chart.set_js` looks the element up by that id when wiring up the D3 chart.
`spectrum_chart.set_js` is the full JavaScript that feeds the chart its spectrum + peak data
— interpolate it raw inside a `<script>` tag, no escaping.

```jinja
<script>{{ assets.D3_JS }}</script>
<script>{{ assets.SpectrumChart_JS }}</script>
<style>{{ assets.SpectrumChart_CSS }}</style>

{% if spectrum_chart.have_spectrum %}
  <div id="{{ spectrum_chart.div_id }}" style="height: 30vw; flex: 1 2; overflow: hidden;" class="SpecChart"></div>
  <script>{{ spectrum_chart.set_js }}</script>
{% endif %}
```

#### Aside: rel-eff vs. spectrum chart asymmetry

The rel-eff (§7.2) and spectrum-chart embeddings look superficially similar but differ in who
owns the element id:

- **Rel-eff**: the *template* picks the id and passes it to `new RelEffPlot(id)`. The C++ side
  emits no `div_id` field for rel-eff because there's nothing to coordinate — the id never
  leaves the template.
- **Spectrum chart**: the bundled `spectrum_chart.set_js` has the id baked in (it was generated
  by `D3SpectrumExport::write_js_for_chart` C++-side and references the id via
  `getElementById`, the resize observer, etc.). The C++ side therefore *must* tell the template
  which id to use; it does so via `spectrum_chart.div_id`.

If you're writing a template, the practical rule is: for rel-eff, pick any id you want and use
it consistently; for spectrum chart, always use `id="{{ spectrum_chart.div_id }}"`.

## 8. Example reports

These are runnable templates that render against the JSON in §3.9. They are intentionally
short; consult `std_rel_eff_summary.tmplt.html` and `std_rel_eff_summary.tmplt.txt` for the
fully-featured defaults. Save each into a directory of your choosing and pass that directory
as `include_dir` (or pass an absolute path as `tmplt`) to `render_template`.

> **A pattern you'll see in all three.** Inja's `set_trim_blocks(true)` eats whitespace and a
> newline after every `{% ... %}` close (see §5.7). The trick to preserving a trailing newline
> on a line that ends with a control block is to append an empty interpolation `{{ "" }}` —
> trim only fires on `{% %}` closes, not on `{{ }}` closes. The trick to preserving a space
> between `{% endif %}` and the next char on the same line is to put the space *inside* the
> if/else, not after `{% endif %}`. Both are illustrated below.

### 8.1 Minimal CSV summary

```jinja
nuclide,curve,rel_activity,rel_activity_uncert,enrichment_pct,detector_counts
{% for curve_acts in relative_activities %}{% for nuc in curve_acts.nuclides %}
{{ nuc.name }},{{ curve_acts.curve_index }},{{ nuc.rel_activity }},{{ nuc.rel_activity_uncertainty }},{% if existsIn(nuc, "enrichment") %}{{ pct(nuc.enrichment, 6) }}{% endif %},{% if existsIn(nuc, "detector_counts") %}{{ nuc.detector_counts }}{% endif %}{{ "" }}
{% endfor %}{% endfor %}
```

The trailing `{{ "" }}` keeps the per-nuclide newline alive — without it, every row would
concatenate into one line because the `{% endif %}` closes would each eat their trailing `\n`.

### 8.2 Markdown one-pager

`## ` and `### ` at column zero collide with Inja's line-statement prefix (§5.7), so this
template uses inline `<h2>` / `<h3>` HTML — Markdown allows raw HTML and renders it as
headings. The `{{ "" }}` after the Status line's `{% endif %}` keeps the line-break before
χ²/dof. Inside table rows, the trailing space of each cell sits *inside* the `{% if %}` so
that the cell separators line up after trim eats whitespace.

```jinja
# {{ options.spectrum_title }}

**Status:** {{ status.fail_reason }}{% if not status.success %} — {{ status.error_message }}{% endif %}{{ "" }}
**χ²/dof:** {{ chi2_str }} / {{ dof }} = {{ chi2_per_dof_str }}
**Solver:** {{ timing.duration_str }} ({{ timing.solve_calls }} fcn evals + {{ timing.cov_calls }} for covariance)

<h2>Relative-Efficiency equations</h2>
{% for curve in rel_eff_curves %}
- **{{ curve.name }}** ({{ curve.rel_eff_eqn_type }}): `{{ curve.equation_text }}`
{% endfor %}

<h2>Activities</h2>
{% for curve_acts in relative_activities %}
<h3>Curve {{ curve_acts.curve_index }}: {{ length(curve_acts.nuclides) }} nuclide(s)</h3>

| Nuclide | Rel. Activity | Mass % | Enrichment |
|---|---|---|---|
{% for nuc in curve_acts.nuclides %}
| {{ nuc.name }} | {{ printCompact(nuc.rel_activity, 6) }} ± {{ printCompact(nuc.rel_activity_uncertainty, 6) }} | {% if existsIn(nuc, "total_mass_fraction") %}{{ pct(nuc.total_mass_fraction, 4) }}% {% endif %}| {% if existsIn(nuc, "enrichment") %}{{ pct(nuc.enrichment, 4) }}% {% endif %}|
{% endfor %}
{% endfor %}

{% if length(warnings) > 0 %}
<h2>Warnings</h2>
{% for w in warnings %}
- {{ w }}
{% endfor %}
{% endif %}
```

### 8.3 Plain-text peak listing

The trailing whitespace for each column lives *inside* the if/else so that trim doesn't eat
the inter-column gap. Free-floating peaks render `--` in the nuclide / yield / fit-rel-eff
columns; nuclide-assigned peaks render the actual values.

```jinja
Energy (keV)  Nuclide   Yield     Amplitude        FitRelEff   Cont. range
------------  --------  --------  ---------------  ----------  -----------------
{% for p in peaks %}
{{ printFixed(p.energy, 2) }}        {% if p.is_floating %}--        {% else %}{{ p.nuclide }}    {% endif %}{% if p.is_floating %}--          {% else %}{{ printCompact(p.yield, 3) }}    {% endif %}{{ printCompact(p.amplitude, 6) }} ± {{ printCompact(p.amplitude_uncertainty_percent, 2) }}%   {% if p.is_floating %}--          {% else %}{{ printCompact(p.fit_rel_eff, 6) }}    {% endif %}{{ printFixed(p.continuum.lower_energy, 1) }}–{{ printFixed(p.continuum.upper_energy, 1) }}
{% endfor %}
```

## 9. Suggested LLM prompt

When using this README to generate a new template, an LLM should be told:

> *You are writing a single Inja template (Inja 3.x) that consumes the JSON described in the
> README of `InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/`. Use only the
> features listed in §4. Do not use any feature listed in §5 ("does not work"). Prefer the
> custom callbacks (`printFixed`, `printCompact`, `pct`, `safe_html`, `existsIn`, `length`) over
> Inja's expression arithmetic. Output only the template body — no explanatory prose.*

Then provide:
- The README (this file).
- The user's description of the desired report (purpose, audience, format, layout).
- Optional: the user's preferred row/column conventions or styling.

## 10. Where the C++ side lives

| File | Role |
|---|---|
| [`InterSpec/RelActAutoReport.h`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec/RelActAutoReport.h) | Public API |
| [`src/RelActAutoReport.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/RelActAutoReport.cpp) | `solution_to_json` + Inja env + render dispatch |
| [`InterSpec/BatchRelActAuto.h`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec/BatchRelActAuto.h) / [`src/BatchRelActAuto.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/BatchRelActAuto.cpp) | Batch driver — builds the multi-file `summary_json` consumed by §11 |
| [`src/BatchInfoLog.cpp`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/src/BatchInfoLog.cpp) | Sibling module for shielding/source and peak-fit reports; `printFixed`/`printCompact` live here |

To regenerate the JSON example for any solution, pass `tmplt = "json"` to `render_template`
(or call `solution_to_json(sol).dump(2)` directly). The dump always contains the heavy JS/CSS
asset blobs — pipe through `jq 'del(.assets, .spectrum_chart.set_js)'` (see §3.9) for a
readable view.

## 11. Multi-file batch-summary templates

The single-file template you've been reading about for the last ten sections consumes the JSON
returned by `RelActAutoReport::solution_to_json()` for **one** solved spectrum. When InterSpec
runs a batch (CLI: `--batch-iso-from-nucs --input-file=… --input-file=…` or
GUI: Batch dialog → "Isotopics" tab with multiple input files) it also produces an aggregate
**summary** JSON that wraps each file's per-file JSON in a `Files[]` array.  The aggregate is
consumed by a different family of templates — multi-file summary templates — selected via the
`"html-summary"` / `"summary"` shorthand or by a custom filename passed to:

| | CLI flag | GUI control |
|---|---|---|
| Default summary template | `--summary-report-template=html-summary` (also: `default`, `summary`) | "Custom Summary" report uploader on the Isotopics tab |

The single bundled multi-file template,
`std_multi_file_summary.tmplt.html`, opens with a one-line "Analyzed N files for 'Isotopics by
Nuclides', with M succeeding and K failing" banner, then a per-file activity / enrichment summary
table, then per-file detail sections (spectrum chart + rel-eff chart + isotopics tables).

### 11.1 The multi-file JSON shape

The aggregate `summary_json` that batch templates consume is built by
`BatchRelActAuto::run_in_files()`. The shape is:

```json
{
  "ExemplarFile": "exemplar.n42",
  "ExemplarSampleNumbers": [1],
  "InputFiles": [ "fg1.n42", "fg2.n42", ... ],
  "NumFiles": 15,
  "NumSucceeded": 13,
  "NumFailed": 2,
  "assets": { "D3_JS": …, "SpectrumChart_JS": …, "SpectrumChart_CSS": …, "RelEffPlot_JS": …, "RelEffPlot_CSS": … },
  "Files": [
    {
      "Filepath":     "/abs/path/to/fg1.n42",
      "Filename":     "fg1.n42",
      "ParentDir":    "/abs/path/to",
      "ResultCode":   "Success",                  // see BatchRelActAuto::ResultCode
      "ResultCodeInt": 0,
      "Success":      true,
      "HasErrorMessage": false,
      "ErrorMessage": "",                          // present only when HasErrorMessage
      "HasWarnings":  false,
      "Warnings":     ["…"],                       // present only when HasWarnings
      "status":       { … },                       // ← every key from §3 (solution_to_json)
      "spectrum_title": "…",
      "rel_eff_curves": [ … ],
      "relative_activities": [ … ],
      "peaks": [ … ],
      "spectrum_chart": { … },
      "rel_eff_chart_json": "[…]"
      // NOTE: no per-file `assets` here -- they're lifted to the top-level `assets`.
    },
    { … },
    …
  ],
  "Warnings": [ "…", "…" ]
}
```

| Top-level key | Type | Notes |
|---|---|---|
| `ExemplarFile` | string | Path of the exemplar N42 (or empty if a `--rel-eff-config-file` was used instead) |
| `ExemplarSampleNumbers` | array of int | Optional: present only when sample-number disambiguation was supplied to the batch driver |
| `InputFiles` | array of string | Verbatim copy of the input-file list passed on the command line / dropped on the GUI |
| `NumFiles`, `NumSucceeded`, `NumFailed` | int | Counts. Pre-computed server-side for the multi-file template's banner |
| `assets` | object | The heavy JS/CSS asset blobs (`D3_JS`, `SpectrumChart_JS`, `SpectrumChart_CSS`, `RelEffPlot_JS`, `RelEffPlot_CSS`) — **lifted out of each per-file solution and stored once at the top level** so the multi-file template can `{{ assets.D3_JS }}` etc. without duplicating hundreds of KB across N files. Mirrors the pattern used by `BatchActivity::fit_activities_in_files`. |
| `Files` | array | One entry per input file. **Each entry's body is the per-file JSON** documented in §3 (`solution_to_json` output) — *minus* the `assets` block that was lifted up to the top level — plus the bookkeeping fields below |
| `Warnings` | array of string | Whole-batch warnings emitted by `BatchRelActAuto::run_in_files` |

#### Per-file bookkeeping fields (added to each `Files[i]`)

| Field | Type | Notes |
|---|---|---|
| `Filename` | string | Basename of the input file (use this for table rows / per-file headings) |
| `Filepath` | string | Absolute path that was given to the batch driver |
| `ParentDir` | string | Directory containing `Filepath` |
| `ResultCode` | string | One of `BatchRelActAuto::ResultCode`'s names: `Success`, `NoExemplar`, `CouldntOpenExemplar`, `ExemplarMissingRelActState`, `CouldntOpenStateOverride`, `CouldntOpenInputFile`, `CouldntOpenBackgroundFile`, `ForegroundSampleNumberUnderSpecified`, `BackgroundSampleNumberUnderSpecified`, `ExemplarUsesPhysModelButNoDrf`, `FwhmMethodNeedsDrfButNoneAvailable`, `SolveFailedToSetup`, `SolveFailedToSolve`, `SolveUserCanceled`, `SolveThrewException`, `UnknownStatus` |
| `ResultCodeInt` | int | Numeric form of `ResultCode` |
| `Success` | bool | Convenience: `(ResultCode == "Success")` |
| `HasErrorMessage` | bool | True iff `ErrorMessage` is non-empty |
| `ErrorMessage` | string | Present only when `HasErrorMessage` is true |
| `HasWarnings` | bool | True iff `Warnings` is non-empty |
| `Warnings` | array of string | Per-file warnings (also folded into the top-level `Warnings`) |

The remaining keys inside each `Files[i]` are exactly those documented in §3
(`solution_to_json` output) **except** for the `assets` block (lifted to top-level). So:
`status`, `spectrum_title`, `rel_eff_curves`, `relative_activities`, `peaks`, `rois`,
`energy_calibration`, `options`, `curve_options`, `foreground`, `background`, `detector`,
`timing`, `warnings_html`, `timestamps`, `compile_timestamp`, `spectrum_chart`,
`rel_eff_chart_json`, `have_multiple_rel_eff`, `live_time_s`, `chi2`, `dof`, `chi2_per_dof`,
`chi2_str`, `chi2_per_dof_str`, `pu_corrections`, `ratios` — any pattern from §7 lifts cleanly
into a multi-file template by prefixing the path with `file.` inside a `{% for file in Files %}`
loop, *with the single exception that asset references stay top-level*.

### 11.2 Things to know when authoring a multi-file template

1. **Shared assets live at the top level.** The C++ side
   (`BatchRelActAuto::run_in_files`) lifts the JS/CSS asset blobs out of each per-file
   `solution_to_json` payload and stores **one** copy in `summary_json["assets"]`, then erases
   the per-file `assets` field before pushing each entry into `Files[]`. So the template
   embeds them at the top of the document via:
   ```jinja
   <script>{{ assets.D3_JS }}</script>
   <style>{{ assets.SpectrumChart_CSS }}</style>
   ```
   **Do not** reach into `file.assets` or `Files.0.assets` — those keys don't exist on the
   multi-file path. (For variable indexing into other arrays, use the `at()` callback — see
   §5.1.)

2. **Per-file IDs must be unique.** Spectrum charts use `Files[i].spectrum_chart.div_id` (the
   C++ side already generates a per-call random suffix like `specchart_<16-hex-chars>` so the
   id is unique even when multiple solutions are embedded on one page). Rel-eff plot IDs are
   template-controlled — generate them with `loop.index`:
   ```jinja
   <div id="releffchart_{{ loop.index }}" class="…"></div>
   <script>
     const chart_{{ loop.index }} = new RelEffPlot("releffchart_{{ loop.index }}");
     chart_{{ loop.index }}.setRelEffData({{ file.rel_eff_chart_json }});
   </script>
   ```

3. **Always branch on `file.Success` before touching any per-file fit fields.** Failed files
   still appear in `Files[]` (so `NumFiles` / the summary banner / top-level
   `Files` length stays right) but their `relative_activities`, `peaks`, etc. will be
   default-empty. The bundled template renders them as a "FAILED — ResultCode: ErrorMessage"
   row in the summary table and a single failure paragraph in the per-file section.

4. **Don't `safe_html(file.assets.*)` or `safe_html(rel_eff_chart_json)`.** Those are JS / CSS /
   JSON literals meant to be interpolated raw inside `<script>` / `<style>` tags. Apply
   `safe_html` only to user-supplied strings (`file.spectrum_title`, `file.ErrorMessage`,
   warnings).

### 11.3 Minimal multi-file CSV summary (example)

A short multi-file template — one nuclide row per file per curve — useful when the user wants a
spreadsheet-friendly batch summary:

```jinja
file,curve,nuclide,rel_activity,rel_activity_uncert,enrichment_pct,result_code
{% for file in Files %}{% if file.Success %}{% for curve_acts in file.relative_activities %}{% for nuc in curve_acts.nuclides %}
{{ file.Filename }},{{ curve_acts.curve_index }},{{ nuc.name }},{{ nuc.rel_activity }},{{ nuc.rel_activity_uncertainty }},{% if existsIn(nuc, "enrichment") %}{{ pct(nuc.enrichment, 6) }}{% endif %},{{ file.ResultCode }}{{ "" }}
{% endfor %}{% endfor %}{% else %}
{{ file.Filename }},,,,,,{{ file.ResultCode }}{{ "" }}
{% endif %}{% endfor %}
```

The trailing `{{ "" }}` after each row guards the per-row newline against `set_trim_blocks`
(see §5.7); the `{% else %}` branch emits one empty-cell row for failed files so every input
appears at least once in the output.

### 11.4 The bundled `std_multi_file_summary.tmplt.html`

The production multi-file template ships at
[`InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/std_multi_file_summary.tmplt.html`](https://raw.githubusercontent.com/sandialabs/InterSpec/refs/heads/master/InterSpec_resources/static_text/IsotopicsByNuclidesReportTmplts/std_multi_file_summary.tmplt.html)
and demonstrates the patterns from §11.2 in a complete report:

- `Files.0.assets.*` pulled once at the top of `<head>`;
- header line `Analyzed {{ NumFiles }} files for "Isotopics by Nuclides", with {{ NumSucceeded
  }} succeeding and {{ NumFailed }} failing`;
- a single combined activity / enrichment table whose rows are grouped by file (one bold
  `file-row` separator above each file's nuclide rows);
- per-file detail sections, each with its own spectrum chart + rel-eff chart (using
  `loop.index` for the rel-eff div id) + Pu / activity / mass-ratio tables — i.e. a compact
  version of `std_rel_eff_summary.tmplt.html` with `file.` prefixes added to every key access.

When asked to author a fresh multi-file template, an LLM should follow §11.1's JSON schema and
the four authoring rules in §11.2; if the user wants a particular sub-block (e.g. just the
summary table, or just the per-file sections) the bundled template is a reasonable starting
point to copy from.
