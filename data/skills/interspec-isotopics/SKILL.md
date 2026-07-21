---
name: interspec-isotopics
description: >-
  Performs isotopic ratio analysis (enrichment, burn-up) for Plutonium and Uranium from
  gamma-ray spectroscopy using InterSpec MCP tools. Covers preset selection, configuration,
  relative efficiency fitting, ROI management, and result evaluation. Use when determining
  U-235 enrichment, Pu isotopic composition, and non-U/Pu isotopic relative activity analysis.
metadata:
  author: InterSpec
  version: "0.1"
---

# Isotopics Analysis Workflow

Determine isotopic composition (mass fractions) of nuclear materials from gamma-ray spectroscopy, usually Plutonium and Uranium, but can be used for any nuclide. For non-U/Pu spectra, can be useful to get relative activity ratios without knowing detector efficiency or source geometry/distance/shielding.

## Method Overview

Isotopics analysis uses a **relative efficiency** approach: the efficiency as a function of energy, and activities of each nuclide, are fit by comparing observed peak areas to predictions from a given set of activities, efficiency parameters, and optionally shielding thicknesses. Once optimized, relative activities yield isotopic mass fractions. Distance and absolute detector efficiency are not required.

## Workflow

### Step 1: Analyze Request

Determine the analysis type:
- **Plutonium**: Pu-239, Pu-240, Pu-241, Am-241, etc. Look for peaks in the ~100-780 keV range.
- **Uranium**: U-235/U-238 enrichment. Look for peaks near 186 keV (U-235), 1001 keV (Pa-234m from U-238).
- **Custom**: Other isotopic ratios. Start from blank configuration.

Check prerequisites:
- `get_spectrum_info` - Verify foreground spectrum is loaded
- `get_peaks` - Check for relevant peaks to confirm material type

### Step 2: Select and Load Configuration

Use `list_isotopics_presets` to see available presets, then `load_isotopics_preset`.

**For Plutonium:**
| Preset | When to use |
|---|---|
| HPGe Pu (120-780 keV).xml | **Default**. Wide energy range, most comprehensive. |
| HPGe Pu (120-420 keV).xml | Low-energy region only. Use when low-energy peaks are clean. |
| HPGe Pu (610-775 keV).xml | High-energy only. Use when Am-241 interference at 59.5 keV is severe. |

**For Uranium:**
| Preset | When to use |
|---|---|
| HPGe U (120-1001 keV).xml | **Default** for U-235/U-238 enrichment analysis. |

**For other sources:**
- `reset_isotopics_config` to start blank
- `modify_isotopics_rois` to add ROIs covering the major gamma lines (consider restricting to >~120 keV if sufficient high-energy lines exist)
- `modify_isotopics_nuclides` to add the nuclides of interest

### Step 3: Validate Configuration

Call `get_isotopics_config` to review:
- Verify expected nuclides are present
- Verify ROIs cover relevant peaks
- If source ages are known, use `modify_isotopics_nuclides` to set them

### Step 4: Check for Interferences

Review with `get_peaks` (filter 'analysis') for non-target nuclides with peaks in the analysis energy ranges.

**Common interferences:**
- **Pu analysis**: I-131, Cs-137, other fission products
- **Mixed samples**: U and Pu both present
- **Background**: K-40, Ra-226, Th-232 (usually outside ROI ranges)

If significant interference is found, add the interfering nuclide to the configuration with `modify_isotopics_nuclides`, or adjust ROIs with `modify_isotopics_rois`.

### Step 5: Execute Calculation

Call `perform_isotopics_calculation` (no parameters needed - uses current configuration).

Returns:
- Mass fractions with uncertainties for each isotope
- Relative activities
- Ages (if fitted)
- Quality metrics (chi-squared per DOF)
- Pu-242 corrected values (for Pu analysis, if available)

### Step 6: Evaluate Results

**Overall chi2/DOF:**
| Value | Assessment |
|---|---|
| 0.7 - 1.5 | Good fit |
| < 5.0 | Acceptable |
| 5.0 - 10.0 | Possibly acceptable (especially high-statistics spectra) |
| > 10.0 | Poor - solution likely not acceptable |

**Per-ROI evaluation:**
- Check `roi_info -> chi2_per_dof` for individual ROIs to find problem regions
- High chi2 in a single ROI may indicate missing source, wrong continuum type, or interference
- For high-statistics ROIs, try changing continuum type (e.g., to stepped continuum)

**Observed efficiencies:**
- Check `observed_efficiencies -> deviation_from_solution_sigma`
- Under-prediction (negative deviation) may indicate missing source or too much shielding
- Over-prediction (positive deviation) may indicate shielding needed, especially if low-energy peaks over-predict and high-energy under-predict

**If results are poor:**
- Try different peak skew model via `modify_isotopics_options`
- Adjust ROI extents, remove minor ROIs, or add ROIs with more significant peaks
- Try different relative efficiency curve settings (`modify_isotopics_curve_settings`)

### Step 7: Report Results

**For U/Pu**, report mass fractions with uncertainties:
```
Uranium Isotopics Results:
  U-235: 20.0 +/- 1.0 %
  U-238: 78.1 +/- 1.5 %
  U-234: 0.8 +/- 0.01 %

Quality: chi2/DOF = 1.2 (good fit)
```

**For non-U/Pu**, report relative activities with uncertainties.

Include:
- Quality metrics (chi2/DOF overall and any problematic ROIs)
- Age information if fitted
- Pu-242 corrections if applicable
- Caveats about interferences or limitations

## Configuration Modification Tools

- `modify_isotopics_nuclides` - Add/remove/update nuclides. Set ages, age fitting bounds, activity constraints.
- `modify_isotopics_rois` - Add/remove/update energy ROIs. Set continuum type, range behavior.
- `modify_isotopics_curve_settings` - Change relative efficiency equation type (LnXLnY, LnX, LnY, FramPhysicalModel), order, self-attenuation, external shielding.
- `modify_isotopics_options` - Global options: energy cal fit mode, FWHM form, peak skew type.
- `modify_isotopics_constraints` - Activity ratio and mass fraction constraints between nuclides.
- `get_isotopics_config_schema` - Full schema with valid values for all configuration fields.

## Glossary

- **ROI**: Energy range used for isotopics analysis (distinct from peak ROI).
- **Mass fraction**: Percentage of total mass from each isotope.
- **Relative activity**: Activity ratio between isotopes (not absolute activity).
- **Age**: Time since purification of U or Pu (for Pu: when Am-241 was separated from Pu-241).
- **Pu-242 correction**: Empirical correction for Pu-242, which has no usable gamma lines.
