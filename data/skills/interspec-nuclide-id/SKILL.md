---
name: interspec-nuclide-id
description: >-
  Identifies radioactive sources (nuclides, x-ray elements, nuclear reactions) in gamma-ray
  spectra using InterSpec MCP tools. Covers systematic peak investigation, source validation
  via energy matching and relative efficiency, the Ultimate Parent rule, and NORM background
  handling. Use when identifying unknown peaks or performing full-spectrum source identification.
metadata:
  author: InterSpec
  version: "0.1"
---

# Nuclide Identification Workflow

Systematic workflow for identifying the source of photopeaks in a gamma-ray spectrum.

## Goal

Identify the source (nuclide, x-ray, nuclear reaction) of every significant peak in the spectrum. Work until all peaks are identified or rigorously determined to be unidentifiable.

## The Ultimate Parent Rule

Do not assign peaks to short-lived daughters. Determine the long-lived parent:
- Pb-214, Bi-214 peaks -> assign to **Ra-226** (or U-238)
- Ac-228, Pb-212, Tl-208 peaks -> assign to **Th-232**
- Pa-234m, Th-234 peaks -> assign to **U-238**
- Ba-137m -> assign to **Cs-137**

Only assign to the daughter if there is evidence of chemical separation.

## Workflow

### Step 1: Assess Spectrum State

Call these tools first:
- `get_loaded_spectra` - What spectra are available?
- `get_spectrum_info` - Detector type, energy range, live time
- `get_identified_sources` - Sources already assigned to analysis peaks
- `get_peaks` with `{"filter": "analysis"}` - Current analysis peak state
- `get_automated_id_results` - Hints from automated algorithms (may be wrong)

Use `source_info` on any already-identified sources to learn about commonly associated nuclides.

**Decision**: Proceed to Step 2 (Get Unidentified Peaks). If a background spectrum is loaded, NORM peaks are handled automatically. If no background, NORM peaks are handled on-demand during investigation when encountered.

### Step 2: Get Unidentified Peaks

Call `get_peaks` with `{"filter": "unidentified", "maxResults": 3}` (or more).

Select the most distinct or easiest-to-identify peak (highest amplitude, highest energy, or tightest FWHM). If all peaks have been investigated, proceed to Step 6a (Check Non-Peak Signatures).

**Important**: After fitting any peaks (`add_analysis_peak` or `add_analysis_peaks_for_source`), the unidentified peaks list becomes stale. Call `get_peaks` with the 'unidentified' filter again.

### Step 3: Investigate Peak

For the selected peak energy:
1. `search_sources_by_energy` with the peak energy. Prefer single-energy searches. The `profile_score` metric (higher = better) accounts for spectrum fit.
2. Check for escape peaks: `escape_peak_check`
3. Check for sum peaks: `sum_peak_check`
4. Look for parent/daughter relationships (e.g., if 351.9 keV found, consider Ra-226).

If a plausible candidate is found, proceed to Step 4. If not, return to Step 2 and skip this peak.

**NORM Detection (no background only):** If no background spectrum is loaded and NORM peaks have not yet been fit, check if the peak matches common NORM energies:
- K-40: 1460.8 keV
- Ra-226 chain: 609.3, 1120.3, 1764.5, 295.2, 351.9 keV
- Th-232 chain: 583.2, 2614.5, 911.2, 968.9, 338.3, 238.6 keV

If the peak matches NORM, or NORM interference is suspected (e.g., 186 keV where Ra-226 could interfere with U-235), fit all NORM peaks at once:

```
add_analysis_peaks_for_source({source: null, options: ["FitNormPeaks"]})
```

Then check NORM ratios (Ra-226 and Th-232 activities within ~1 order of magnitude, 609.3 keV area > 583.2 keV, K-40 present at 1460.8 keV) and return to Step 2.

### Step 4: Validate Candidate

Before accepting a source identification:

1. **Context**: Call `source_info` on the candidate. Check for associated nuclides, common mis-IDs, and whether it makes sense given other identified sources.

2. **Energy fit**: Source energy must be within ~1 FWHM of the peak mean.

3. **Corroborating peaks**: Use `source_photons({source: "...", sortOrder: "intensity", maxResults: 25})` to check if other significant gammas from this source have corresponding peaks in the spectrum. Consider:
   - Attenuation may suppress low-energy lines
   - Low branching ratio lines may be below detection threshold
   - Nearby stronger peaks may obscure weaker lines
   - The `prominent_energies_keV` field of the `source_photons` result gives the most characteristic energies for a quick check

4. **Ultimate Parent**: If the candidate is a short-lived daughter, determine the parent.

**Example**: Searching on a peak at 186 keV, top candidates might be Ra-226 and U-235. Check for U-235 corroborating peaks at 205.3 and 163.4 keV. If present, U-235 is confirmed. Note that Ra-226 may also be present; its 295.2 and 351.9 keV peaks should be much larger than 186 keV if Ra-226 is the source (higher branching ratios).

If validation passes, proceed to Step 5. If it fails, try another candidate or return to Step 2.

### Step 5: Add Source Peaks

1. **Create checkpoint first**: `create_peak_checkpoint({name_hint: "before_Co60"})` so you can roll back.
2. **Add peaks**: Use `add_analysis_peaks_for_source` (preferred) to fit all detectable peaks for the source at once. This produces better fits than adding peaks individually.
   - **Default** (no options): If the source already has analysis peaks, they are replaced. New ROIs are trimmed to avoid existing ROIs from other sources.
   - **If a source gamma overlaps an existing ROI from a different source** (e.g., U-235 186 keV overlapping Ra-226 186 keV): use `options: ["RefitInterferingAnalysisPeaks"]` so both sources are properly fit in the shared region.
   - **If you want to leave all existing peaks untouched**: use `options: ["DoNotUseExistingRois"]` — source gammas overlapping existing ROIs will be skipped.
3. **Verify**: Check that peaks were added with correct source assignment.

If the source has >2 peaks, validate with relative efficiency (Step 5a).

#### Step 5a: Check Relative Efficiency

Call `peak_based_relative_efficiency` with a single nuclide (e.g., `sources: ["Co60"]`).

Setup:
- Set `eqn_order` to at most (number_of_peaks - 1). E.g., 3 peaks -> order 2.
- If a peak near 511 keV is not among the top ~3 gamma lines for the source, exclude it via `exclude_peak_energies`.
- Set `background_subtract: true` if source peaks overlap with background peaks.

Interpretation:
- **Good**: All `residual_sigma` < 5.
- **Acceptable** if `residual_sigma` > 5 but: peak is small next to a large one, near 511 keV, low energy vs high energy lines (shielding effects), or secular equilibrium issues.
- **Poor**: Try different `eqn_form` values (LnX, LnY, LnXLnY, FramEmpirical). Try excluding problematic small peaks.
- **Reject**: If amplitudes are physically impossible, use `restore_peaks_to_checkpoint` to revert and reconsider.

### Step 5b: Check Associated Sources

After adding a source, use `source_info` to check for associated nuclides:
- Does the source imply a parent or chain member? (e.g., Ra-226 -> check for U-235; Cs-137 -> check for Cs-134 indicating reactor origin)
- If an associated nuclide is plausible and has detectable peaks, validate it (back to Step 4).

Then return to Step 2 for the next unidentified peak.

### Step 6: Check Non-Peak Signatures

Before finalizing, check for radiation signatures that do not manifest as discrete photopeaks:

**Neutron checks:**
1. Check `get_spectrum_info` for `neutronCPS` field. If background is loaded, check `neutronCPS_excess_sigma` (>2 indicates significant neutron excess).
2. Check continuum above 2614.5 keV using `get_counts_in_energy_range` (~2640-3000 keV) for neutron-induced elevation. Caveat: less diagnostic if Co-56 or similar high-energy emitters are identified.
3. Look for neutron reaction peaks among identified sources: H(n,g) 2223.2 keV, Fe(n,n') 846.8 keV, B-10(n,a) 477.6 keV, Ge(n,n') inelastic scatter.
4. If 477 keV region shows excess, use `get_spectrum_image` (430-530 keV) to check for Doppler-broadened triangular B-10 feature.
5. If spectrum extends above 4 MeV, check 4.2-4.7 MeV range for alpha-Be (AmBe/PuBe) signature using `get_counts_in_energy_range` and `get_spectrum_image`.
6. If detector is HPGe, check low-energy region (~50-200 keV) with `get_spectrum_image` for "ski-slope" continuum indicating neutron interactions in the germanium crystal.

**Alpha checks:**
- 4.44 MeV broad feature indicates alpha-Be neutron source (AmBe/PuBe).
- Note any recoil-broadened peaks observed during investigation.

**Beta source check:**
- If foreground is clearly elevated over background but relatively few peaks were identified (e.g., only 511 keV and/or small NORM peaks relative to total counts), use `get_spectrum_image` below ~2 MeV to check if a smooth exponential continuum dominates. This indicates beta sources (Sr-90/Y-90, P-32, Ru-106/Rh-106).
- Check `get_counts_in_energy_range` in broad bands (100-500, 500-1000, 1000-2000 keV) comparing foreground to background — broad elevation without peaks suggests beta source.

Note findings for the final summary, then proceed to Step 7.

### Step 7: Finalize

Summarize findings:
- List all identified sources (ultimate parents) with confidence levels
- List peaks that remain unidentified, with explanation of why they could not be resolved
- Note any non-peak signatures detected: neutron excess, neutron reaction indicators, alpha-Be features, HPGe ski-slope, beta continuum, Doppler-broadened features
- Note any anomalies (unusual NORM ratios, unexpected sources, failed fits)

Do not add unidentified peaks with source "Unknown" to the analysis state.
