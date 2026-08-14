---
name: interspec-mcp
description: >-
  Gamma-ray spectrum analysis via the InterSpec MCP server. Provides tools for loading spectra,
  fitting peaks, identifying radioactive sources (nuclides, x-rays, reactions), calculating
  decay, determining detection limits (MDA), fitting source activities and shielding, performing
  isotopics (enrichment/burn-up) analysis, and managing energy calibration. Use when analyzing
  radiation spectra, identifying nuclides, performing nuclear calculations, or working with
  InterSpec.
metadata:
  author: InterSpec
  version: "0.1"
---

# InterSpec MCP: Gamma-Ray Spectrum Analysis

InterSpec is a nuclear radiation spectral analysis application. The MCP server exposes tools for programmatic spectrum analysis, nuclide identification, activity fitting, isotopics, and nuclear data queries.

## Core Principles

1. **Tool-first**: Always use tools to access spectrum data and nuclear information. Tool calls are cheap; assumptions are expensive.
2. **Expert audience**: Assume the user is a knowledgeable spectroscopist. Focus on data-driven conclusions.
3. **Step-by-step reasoning**: Articulate your analytical process. Explain why you take certain steps.
4. **Terminology precision**: Distinguish between **detected peaks** (from automated search) and **analysis peaks** (in the internal state used for fitting). Distinguish **CPS** (counts per second, for comparing spectra) from raw counts.

## Key Concepts

- **Analysis peaks**: Peaks in the application's internal state used for nuclide fitting and activity analysis. Must have an attributed source to be useful. Managed via `add_analysis_peak`, `add_analysis_peaks_for_source`, `edit_analysis_peak`.
- **Detected peaks**: Peaks found by automated peak search. May be incomplete; additional peaks may be fittable.
- **ROI (Region of Interest)**: One or more peaks sharing a continuum. Defined by energy bounds and a continuum type.
- **CPS**: Counts per second. Use CPS (not amplitude) when comparing peak intensities between spectra with different live times.
- **FWHM**: Full width at half maximum (~2.355 * sigma for a Gaussian). Use `get_expected_fwhm` to determine expected width.
- **Source types**: Nuclides (e.g., Co60, Cs137), fluorescence x-ray elements (e.g., Pb, U), nuclear reactions (e.g., H(n,g)).
- **Ultimate Parent rule**: Assign peaks to the long-lived parent, not short-lived daughters. E.g., assign Pb-214/Bi-214 peaks to Ra-226, not to the daughter directly. See the [interspec-nuclide-id](../interspec-nuclide-id/SKILL.md) sub-skill.
- **NORM**: Naturally Occurring Radioactive Material (K-40, Ra-226/U-238 chain, Th-232 chain). Common background in most spectra.
- **DRF**: Detector Response Function. Describes detector efficiency as a function of energy. Required for activity fitting and MDA calculations.

## Tool Reference

### Spectrum Management
- `get_loaded_spectra` - List currently loaded spectrum types (Foreground, Background, Secondary)
- `get_spectrum_info` - Detailed metadata: live time, energy range, calibration, detector type. Pass `fileInfo: true` for full file details.
- `load_spectrum_file` - Load a spectrum file from the filesystem. Load foreground first if loading a pair (background may be cleared when foreground changes).
- `get_counts_in_energy_range` - Photon counts in an energy range with statistical foreground/background comparison.

### Peak Detection and Management
- `get_peaks` - Peaks in a spectrum, grouped by ROI. `filter` selects the subset: `all` (default, automated search + analysis peaks), `analysis` (internal analysis state - user-selected or agent-added), `unidentified` (foreground peaks without source assignments, elevated over background, sorted by amplitude), or `elevated_above_background`.
- `get_identified_sources` - Unique sources currently assigned to analysis peaks.
- `add_analysis_peak` - Fit and add a single peak at a specified energy, optionally assigning a source.
- `add_analysis_peaks_for_source` - Fit all detectable peaks for one or more sources at once. Much better than fitting individually for multi-line sources. Use `options: ["FitNormPeaks"]` with `source: null` for NORM background. By default, existing peaks for the same source are replaced and new ROIs are trimmed to avoid other sources' ROIs. Use `options: ["RefitInterferingAnalysisPeaks"]` when a source gamma overlaps an existing ROI from a different source (e.g., U-235/Ra-226 at 186 keV). Use `options: ["DoNotUseExistingRois"]` to skip source gammas that overlap any existing ROI.
- `edit_analysis_peak` - Everyday peak edits in one call (centroid, fwhm, amplitude, fit-for flags, ROI bounds, source, color, label, usage flags), with automatic refit. Also supports structural actions (delete, split, merge).
- `set_peak_shape` - Expert peak-shape surgery: skew type/params, continuum type/coefficients, uncertainty overrides.
- `get_expected_fwhm` - Expected peak width at a given energy.

### Source Identification
- `search_sources_by_energy` - Search for nuclides/x-rays/reactions matching specified energies, ranked against the loaded spectrum. `profile_score` (higher = better) is the best match indicator. Prefer single-energy searches.
- `sources_with_gammas_near_energy` - Catalog lookup: common field sources with a characteristic gamma near an energy (optional `window`); does not consider the loaded spectrum.

### Nuclear Data
- `source_info` - Half-life, decay modes, analyst notes, associated sources, common mis-IDs, source categories.
- `source_photons` - Photon lines (energy/intensity table) for a source, plus `prominent_energies_keV` (its most characteristic energies). Supports sorting by energy or intensity, energy range filtering, cascade sums.
- `nuclide_decay_chain` - Full decay chain with branching ratios.
- `decay_calculator` - Calculate radioactive decay over time. Supports forward and back-decay, date ranges, progeny inclusion.
- `get_automated_id_results` - Results from on-board and GADRAS automated ID algorithms (may be wrong/incomplete).

### Peak Validation
- `escape_peak_check` - Check if a peak is a single/double escape peak of a higher-energy parent.
- `sum_peak_check` - Check if a peak is a cascade-sum or random-sum (pile-up) peak.
- `beta_continuum_check` - Deterministic check for a pure-beta-emitter (bremsstrahlung) continuum; background-subtracts the foreground and grades the net continuum shape (requires a background loaded). Returns candidate pure-beta nuclides consistent with the observed termination energy.
- `peak_based_relative_efficiency` - Validate source ID by checking if peak areas are consistent with a smooth relative efficiency curve. Run on one nuclide at a time. See [interspec-nuclide-id](../interspec-nuclide-id/SKILL.md).

### Detection and Dose
- `currie_mda_calc` - Minimum Detectable Activity using Currie/ISO 11929 methodology. Also checks if a peak appears present.
- `photopeak_detection_efficiency` - Fraction of photons that contribute to a peak, accounting for shielding, distance, and detector efficiency.
- `calculate_dose` - Gamma dose rate from a nuclide at a distance (ANSI/ANS-6.1.1-1991, AP orientation).

### Detector Efficiency
- `available_detector_efficiency_functions` - List available DRFs.
- `load_detector_efficiency_function` - Load a DRF by name, path, formula, or energy/efficiency pairs.
- `detector_efficiency_function_info` - Info about the currently loaded or a named DRF.

### Materials and Shielding
- `get_materials` - List available shielding materials.
- `get_material_info` - Material properties (density, atomic number, composition).

### Activity and Shielding Fitting
- `activity_fit` - Fit activities using current GUI configuration (`from_app_state` mode).
- `activity_fit_one_off` - Custom fit with full specification (distance, geometry, shielding, peaks). Does not modify GUI state.
- `get_shielding_source_config` - Current fit configuration from app state.
- `modify_shielding_source_config` - Modify fit configuration (add/remove shielding, set distance/geometry, source ages).
- `mark_peaks_for_activity_fit` - Batch mark/unmark peaks for use in activity fitting.
See the [interspec-activity-fit](../interspec-activity-fit/SKILL.md) sub-skill for detailed workflow.

### Isotopics Analysis
- `list_isotopics_presets` - Available preset configurations (Pu, U variants).
- `load_isotopics_preset` - Load a preset configuration.
- `get_isotopics_config` / `reset_isotopics_config` - Get or clear current configuration.
- `perform_isotopics_calculation` - Execute isotopics analysis using current configuration.
- `modify_isotopics_nuclides` / `modify_isotopics_rois` / `modify_isotopics_curve_settings` / `modify_isotopics_options` / `modify_isotopics_constraints` - Modify configuration details.
- `get_isotopics_config_schema` - Full schema and valid values for all fields.
See the [interspec-isotopics](../interspec-isotopics/SKILL.md) sub-skill for detailed workflow.

### Energy Calibration
- `fit_energy_calibration` - Fit energy calibration using known peaks. Supports gain-only fit with propagation.
- `save_energy_cal_checkpoint` / `restore_energy_cal_checkpoint` - Checkpoint and restore energy calibration state.
See the [interspec-energy-cal](../interspec-energy-cal/SKILL.md) sub-skill for detailed workflow.

### Checkpoints
- `create_peak_checkpoint` - Snapshot current analysis peaks for undo capability.
- `restore_peaks_to_checkpoint` - Restore peaks to a previous checkpoint.

### Relative Efficiency
- `peak_based_relative_efficiency` - Validate source IDs by checking peak area consistency.
- `get_rel_act_manual_state` - Current relative efficiency configuration state.

## Quick-Start Workflows

### Identify Sources in a Spectrum
1. Check what's loaded (`get_loaded_spectra`, `get_spectrum_info`).
2. If background is present, use `get_peaks` with `{"filter": "unidentified"}` to find elevated foreground peaks.
3. If no background, fit NORM first: `add_analysis_peaks_for_source({source: null, options: ["FitNormPeaks"]})`.
4. For each unidentified peak, use `search_sources_by_energy`, validate with `source_info` and corroborating peaks, then `add_analysis_peaks_for_source`.
5. See the [interspec-nuclide-id](../interspec-nuclide-id/SKILL.md) sub-skill for the full workflow.

### Fit Source Activity
1. Ensure peaks are assigned to nuclides and a DRF is loaded (`detector_efficiency_function_info`).
2. Use `activity_fit_one_off` for a quick custom fit, or configure and run `activity_fit` via app state.
3. See the [interspec-activity-fit](../interspec-activity-fit/SKILL.md) sub-skill for details.

### Determine Isotopics (Enrichment/Burn-up)
1. Load appropriate preset (`list_isotopics_presets`, `load_isotopics_preset`).
2. Validate configuration and run `perform_isotopics_calculation`.
3. See the [interspec-isotopics](../interspec-isotopics/SKILL.md) sub-skill.

### Adjust Energy Calibration
1. Create checkpoints (peaks + energy cal).
2. Find NORM peaks (K-40 at 1460.8 keV, Th-232 at 2614.5 keV).
3. Run `fit_energy_calibration`, evaluate, restore peak checkpoint.
4. See the [interspec-energy-cal](../interspec-energy-cal/SKILL.md) sub-skill.

### Calculate Decay
Use `decay_calculator` with nuclide, activity, and time duration or date range. Use `include_progeny: true` to see daughter activities.

### Determine Detection Limits
Use `currie_mda_calc` with the energy of interest. Add nuclide and distance for activity-based MDA.

## Common Pitfalls

- **Comparing spectra**: Always compare CPS, not raw counts or amplitude, between foreground and background.
- **Multi-line sources**: Use `add_analysis_peaks_for_source` instead of fitting peaks individually - it produces much better fits for sources with multiple gamma lines.
- **Short-lived daughters**: Do not assign peaks to short-lived daughters (Pb-214, Bi-214, Ac-228). Assign to the ultimate parent (Ra-226, Th-232) unless chemical separation is known.
- **Background peaks**: You can add background peaks to analysis state, but do not fit their activity as if they were source peaks.
- **Detector efficiency**: Activity fitting and MDA calculations require a loaded DRF. Check with `detector_efficiency_function_info` first.
- **Energy calibration**: If peaks appear systematically shifted from expected energies, calibrate first before nuclide ID.
- **Peak search completeness**: The automated peak search may miss some peaks. Additional peaks may be fittable in the data.

For detailed domain knowledge (NORM ratios, background triage methodology, validation criteria), see [references/domain-knowledge.md](references/domain-knowledge.md).
