---
name: interspec-energy-cal
description: >-
  Adjusts energy calibration of gamma-ray spectra using NORM peaks (K-40, Th-232) via
  InterSpec MCP tools. Covers checkpoint management, calibration peak selection and
  validation, gain-only fitting, and multi-spectrum calibration. Use when peaks appear
  shifted from expected energies or when the user requests energy calibration adjustment.
metadata:
  author: InterSpec
  version: "0.1"
---

# Energy Calibration Workflow

Adjust the energy calibration of loaded spectra so peak positions match their known photopeak energies. Uses NORM peaks as calibration anchors.

## When to Use

- Peaks appear systematically shifted from expected energies
- User explicitly requests calibration
- Before nuclide ID, to start with better peak energies

## Calibration Strategy

1. **Background first**: If a background spectrum is loaded, calibrate it first and propagate to foreground/secondary.
2. **Then each spectrum independently**: Even after propagation, calibrate the foreground (and secondary if loaded) individually. Each spectrum has its own small calibration differences.
3. **NORM peaks**: Use K-40 (1460.83 keV) and/or Th-232 (2614.51 keV) as anchors, unless a spectrum contains a known source.
4. **Gain-only fit**: Prefer using a single peak for a gain-only (1-parameter) calibration fit - although if multiple peaks from a known source (ex Eu152) are known, then the gain and offset may be fit using many peaks if really desired.

## Assumptions

- Current calibration is roughly correct but may be off by up to ~15 keV
- Only polynomial and full-range-fraction calibrations can be fit (LowerChannelEdge is not supported)

## Workflow

### Step 1: Assess Spectra

- `loaded_spectra` - What spectra are available?
- `get_spectrum_info` for each - Check calibration type and energy range

### Step 2: Create Checkpoints

**Always create both checkpoints before any changes:**
- `create_peak_checkpoint({name_hint: "pre_energy_cal"})` - Snapshot analysis peaks
- `save_energy_cal_checkpoint({name_hint: "pre_energy_cal"})` - Snapshot energy calibration

### Step 3: Find Calibration Peaks

Start with background spectrum if present, otherwise foreground.

1. Check existing analysis peaks (`get_analysis_peaks`) for peaks near 1460.83 keV (K-40) or 2614.51 keV (Th-232) with correct source assignment.
2. If no suitable analysis peak exists, use `get_detected_peaks` to find candidates. Search within ~15 keV of expected energy (calibration may be off).
3. If using a NORM peak, pick **one peak** for gain-only fit:
   - Prefer Th-232 at 2614.51 keV if `numSigma > 12` (higher energy = better gain constraint)
   - Otherwise use whichever peak (K-40 or Th-232) has larger `numSigma`
   - If user specified a different peak (e.g., Cs-137 at 661.7 keV), use that instead
4. If the chosen peak is not already an analysis peak, add it: `add_analysis_peak({energy: <peak_energy>, source: "K40", specType: "<current_spectrum>"})`. Verify the correct gamma was assigned by checking `source.photonEnergy` in the result.

### Step 4: Validate Peak Identity

If the detected peak is more than ~2 keV from the expected NORM energy, validate before using:

Use `search_sources_by_energy` at the peak mean energy:
- Near 1460 keV: rule out Eu-152 at 1457.64 keV
- Near 2614 keV: rule out Co-56 at 2598.5 keV
- Check if candidates have other lines that do not correspond to detected peaks

Only proceed if confident the peak is from K-40 or Th-232. If ambiguous, try the other NORM peak or report failure.

**Note for Th-232**: The 2614 keV peak comes from Tl-208 decay to Pb-208. Any Th-232 progeny between Th-232 and Tl-208 can be ignored as alternative sources.

### Step 5: Fit Calibration

Call `fit_energy_calibration` with the validated peak energy.

**For background spectrum**: Use `propagate_to` option of the tool-call to propagate the calibration change to foreground and secondary (if loaded). Background and foreground often share calibration offsets.

**For foreground/secondary** (calibrated individually after background): Do NOT use `propagate_to`. Changes apply only to that spectrum.

### Step 6: Evaluate Result

Examine the per-peak results from `fit_energy_calibration`:
- If `sum_abs_diff_after_keV < sum_abs_diff_before_keV`, calibration improved
- If peak difference after fit is < ~0.5 keV or less than 0.1 FWHM, calibration is good
- If calibration made things worse, restore using `restore_energy_cal_checkpoint`

### Step 7: Calibrate Next Spectrum

**Mandatory**: You must calibrate each loaded spectrum independently, even if the background calibration was propagated. Propagation gives an approximate correction; each spectrum needs its own fit.

Order: Background -> Foreground -> Secondary (if loaded).

For each remaining spectrum, repeat Steps 3-6.

### Step 8: Finish

When all spectra have been independently calibrated:

1. **Restore peak checkpoint**: `restore_peaks_to_checkpoint` with the checkpoint name from Step 2. This removes any peaks added during calibration so the agent does not permanently alter the analysis state.
2. **Do NOT restore the energy cal checkpoint** - the calibration changes should be kept.
3. Provide a summary of calibration changes (before/after peak energy differences for each spectrum).

### If Calibration Fails

If no suitable NORM peaks are found, peak identity is ambiguous, or calibration made things worse:

1. `restore_peaks_to_checkpoint` - Undo peak changes
2. `restore_energy_cal_checkpoint` - Undo calibration changes
3. Report why calibration failed

## Peak Selection Summary

| Condition | Peak to use |
|---|---|
| Th-232 2614 keV has numSigma > 12 | Th-232 2614.51 keV |
| Both present, Th-232 numSigma <= 12 | Whichever has larger numSigma |
| Only K-40 visible | K-40 1460.83 keV |
| User specifies a peak/source | Use that peak |
| No NORM peaks found | Report failure |
