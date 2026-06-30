---
name: interspec-activity-fit
description: >-
  Determines radioactive source activities and shielding parameters from identified peaks
  using InterSpec MCP tools. Covers fitting modes, geometry configuration, shielding models,
  trace sources, age fitting, and result interpretation. Use when fitting source activity,
  determining shielding thickness, or estimating dose-relevant quantities.
metadata:
  author: InterSpec
  version: "0.1"
---

# Activity and Shielding Fitting Workflow

Determine the activity (Bq or Ci) and shielding parameters for identified radioactive sources.

## Prerequisites

Before fitting, verify:
1. **Detector efficiency loaded**: Call `detector_efficiency_function_info` (no args). If not loaded, use `available_detector_efficiency_functions` then `load_detector_efficiency_function`.
2. **Peaks assigned to nuclides**: Call `get_analysis_peaks` and verify peaks have source assignments.
3. **Peaks marked for fitting**: Check the `useForShieldingSourceFit` attribute on analysis peaks. Use `mark_peaks_for_activity_fit` or `edit_analysis_peak` with `useForShieldingSourceFit` to adjust.

## Detector Efficiency Types

The loaded DRF type fundamentally affects how fits are configured. Use `detector_efficiency_function_info` to check.

### Far-Field (FarField) DRFs
Characterize the detector's intrinsic or absolute efficiency as a function of energy, independent of source geometry. The source-to-detector distance and detector face diameter are used to compute the geometric solid angle.

- **FarFieldIntrinsic**: Efficiency is the probability that a photon hitting the detector face produces a full-energy count. Requires `distance` and detector diameter (stored in the DRF) to compute activity.
- **FarFieldAbsolute**: Efficiency is measured at a specific characterization distance. Requires `distance` to scale.

With far-field DRFs, you **must specify distance** in the fit configuration. Fitted activity is in absolute units (Bq, Ci).  These DRFs are valid to be used at any distance, although very close up (less than ~10 cm), the far-field approximation breaks down, so accuracy is reduced.

### Fixed-Geometry (FixedGeom) DRFs
Characterize the detector for a specific, fixed measurement geometry (e.g., a sample sitting directly on the detector face in a specific container). The efficiency already accounts for the geometry, so no distance parameter is needed.

- **FixedGeomTotalAct**: Activity result is total activity in the sample.
- **FixedGeomActPerCm2**: Activity result is per unit area (Bq/cm2).
- **FixedGeomActPerM2**: Activity result is per unit area (Bq/m2).
- **FixedGeomActPerGram**: Activity result is per unit mass (Bq/g).

With fixed-geometry DRFs, **do not specify distance** (it is not meaningful). Shielding and geometry options are generally not applicable. The fit yields activity directly in the units implied by the DRF type.

### Impact on Fit Configuration
- **Far-field**: Specify distance, geometry, and shielding as needed. This is the most common case.
- **Fixed-geometry**: Omit distance. Use `activity_fit_one_off` with just `peak_energies` and `geometry: "Spherical"` (or use `activity_fit` from app state). The result activity units are determined by the DRF type.

## Fitting Modes

### `activity_fit` (from app state)
Uses the GUI configuration. Best when:
- User has set up the fit interactively
- You want results reflected in the GUI
- Iterating on an existing configuration

### `activity_fit_one_off` (custom)
Fully specified fit, does not modify GUI state. Best for:
- Quick exploratory fits
- "What-if" scenarios
- When you want full control over all parameters

## Configuration

### Default Assumptions
- **Geometry**: Spherical (point source) unless extended source is mentioned
- **Distance**: 100 cm (but return an error if not specified for non-fixed-geometry detectors)
- **Shielding**: None unless user mentions shielding
- **Age fitting**: OFF unless user mentions age or nuclide has rich progeny
- **Trace source activity type**: `per_cm3` unless specified otherwise

### Geometry Selection

| User describes... | Geometry |
|---|---|
| Point source, nothing specific | Spherical |
| Pipe, cylinder, rod | CylinderSideOn or CylinderEndOn |
| Box, plate, rectangular | Rectangular |
| Ground/surface contamination | CylinderEndOn (~1 cm Air shielding with trace source, large radius ~100m) |

### Shielding Configuration

Shielding layers are ordered from innermost (closest to source) to outermost:

```json
{
  "shielding": [
    {"material": "Fe", "radial_thickness": "0.5 cm"},
    {"material": "Pb", "fit_thickness": true}
  ]
}
```

- If `radial_thickness` is specified, it's fixed. If omitted, thickness is fit.
- Use `fit_thickness: true` to explicitly fit even when a starting value is provided.
- Material can be an element symbol (Fe, Pb) or a material name from `get_materials`.

### Trace Sources

For volumetric contamination distributed within a shielding layer:

```json
{
  "shielding": [{
    "material": "Al",
    "radial_thickness": "5 mm",
    "trace_sources": [{
      "nuclide": "Cs137",
      "activity": "1 Bq/cm3",
      "fit_activity": true
    }]
  }]
}
```

Activity types auto-detected from units: `Bq/cm3` -> per_cm3, `Bq/g` -> per_gram, `uCi` -> total_activity.

### Surface Contamination

Model surface or ground contamination as a thin disk of Air with a trace source, using CylinderEndOn geometry:

```json
{
  "peak_energies": [661.7],
  "distance": "100 cm",
  "geometry": "CylinderEndOn",
  "shielding": [{
    "material": "Air",
    "length": "1 cm",
    "radius": "100 m",
    "trace_sources": [{
      "nuclide": "Cs137",
      "activity": "1 Bq/cm3",
      "fit_activity": true
    }]
  }]
}
```

The thin length (1 cm) makes it effectively a surface source. The large radius (100 m) models a wide-area contamination. Adjust radius to match the actual contaminated area (if known), or make large enough that the entire contributing area is effectively captured.

## Peak Selection

### For `activity_fit_one_off`
Specify peak energies directly via `peak_energies`. Only include peaks of the isotopes of interest.

### For `activity_fit` (app state)
Uses peaks with the `useForShieldingSourceFit` flag. To manage:
- `mark_peaks_for_activity_fit` - Batch mark/unmark by energy
- `edit_analysis_peak` with `useForShieldingSourceFit` - Individual peaks

**Important**: Only include peaks of nuclides whose activity you want to fit. Unmark NORM peaks (K-40, Ra-226, Th-232) if fitting non-NORM sources.

## Age Fitting

Enable age fitting when:
- User explicitly requests it
- Nuclide has >=2 progeny peaks (check `num_progeny_peaks` in `get_shielding_source_config`)
- The nuclide's age is potentially fittable (check `age_potentially_fittable` - not allowed for equilibrium-limited nuclides like Cs-137)

Isotopes of the same element typically share age (e.g., U-235 and U-238). Use `same_age_isotopes: true`.

## Executing the Fit

### One-off Example

```json
{
  "peak_energies": [1173.2, 1332.5],
  "distance": "100 cm",
  "geometry": "Spherical",
  "shielding": [
    {"material": "Fe", "radial_thickness": "0.5 cm", "fit_thickness": true}
  ]
}
```

### App State Workflow

1. `get_shielding_source_config` - Review current configuration
2. `modify_shielding_source_config` - Adjust as needed (add shielding, set distance, etc.)
3. `mark_peaks_for_activity_fit` - Ensure correct peaks are marked
4. `activity_fit` - Execute the fit

## Interpreting Results

| Metric | Good | Acceptable | Poor |
|---|---|---|---|
| Chi2/DOF | ~1.0 | 0.5 - 2.0 | > 3.0 |
| Activity uncertainty | < 50% of value | < 100% | > 100% (weak constraint) |

Check per-peak observed vs expected counts. If low-energy peaks are under-predicted, shielding may be needed.

## Report to User

Provide:
- Activity with uncertainty (e.g., "10.5 +/- 0.3 uCi" or "388 +/- 11 kBq")
- Age with uncertainty if fitted (e.g., "20 +/- 3 years")
- Shielding parameters if fitted (e.g., "0.5 +/- 0.1 cm Fe")
- Fit quality (chi2/DOF, number of peaks used)
- Any warnings or recommendations

## Common Scenarios

- **Simple point source**: `geometry: Spherical`, no shielding
- **Shielded source**: Add shielding layers, fit thickness if unknown
- **Contaminated pipe**: `CylinderSideOn` geometry with trace source in material
- **Ground contamination**: `CylinderEndOn`, ~1 cm Air layer with trace source, radius ~100 m
- **Multiple nuclides**: The fit supports multiple nuclides simultaneously

## Error Handling

**Hard errors** (fit won't run):
- No detector efficiency function loaded
- No distance specified for non-fixed-geometry detector
- Peak not found or no nuclide assigned

**Soft warnings** (fit runs with caveats):
- Large uncertainties indicate weak constraints
- Poor chi2/DOF suggests model mismatch
- Start simple (no shielding), add complexity only if needed
