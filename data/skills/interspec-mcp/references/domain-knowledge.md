# Domain Knowledge for Gamma-Ray Spectrum Analysis

## Background Triage

When a background spectrum is loaded, differentiate foreground peaks from background:

- For each detected foreground peak, check for a corresponding background peak at similar energy.
- **Compare CPS, not amplitude.** CPS normalizes for different live times.
  - **Explained by background**: Foreground CPS within ~20% of background CPS. Generally ignorable.
  - **Elevated**: Foreground CPS > ~20% above background CPS. Requires investigation.
  - **Unexplained**: No corresponding background peak. Requires investigation.
- If most background-associated peaks are consistently elevated or suppressed in the foreground, the background spectrum may not be representative.

The `get_peaks` tool with `{"filter": "unidentified"}` handles this automatically, returning only elevated or unexplained foreground peaks.

## NORM Background

Naturally Occurring Radioactive Material present in virtually all spectra:

### Common NORM Nuclides
- **K-40**: Single gamma at 1460.83 keV. Always present.
- **Ra-226** (U-238 chain): Key peaks from daughters Pb-214 (295.2, 351.9 keV) and Bi-214 (609.3, 1120.3, 1764.5 keV).
- **Th-232**: Key peaks from daughters Ac-228 (338.3, 911.2, 968.9 keV), Tl-208 (583.2, 2614.5 keV), Bi-212 (727.3 keV).
- **U-238**: Pa-234m at 1001 keV (weak but distinctive).

### Expected NORM Ratios
- The 609.3 keV Ra-226 (Bi-214) peak area should generally exceed the 583.2 keV Th-232 (Tl-208) peak area, since Ra-226 is typically more abundant and Tl-208 has only ~30% branching from the Th-232 chain.
- K-40 at 1460.8 keV should be present with amplitude reasonable relative to thorium/radium lines.
- Ra-226 and Th-232 activities should be comparable (within roughly an order of magnitude) for typical environmental NORM.
- U-238 activity should be consistent with Ra-226 if in secular equilibrium (disequilibrium is common in some materials).

### Fitting NORM Peaks
Use `add_analysis_peaks_for_source` with `{source: null, options: ["FitNormPeaks"]}` to fit all common NORM peaks at once. This uses K-40, Ra-226 (age ~19 days), Th-232 (age ~28.75 years), and U-238 (age ~120.5 days).

## The Ultimate Parent Rule

Do not assign peaks to short-lived daughters. Determine the ultimate parent:

| If you see peaks from... | The source is likely... | Unless... |
|---|---|---|
| Pb-214, Bi-214 | Ra-226 or U-238 | Chemically separated (e.g., Ra-226 sources) |
| Ac-228, Pb-212, Tl-208, Bi-212 | Th-232 | Chemically separated |
| Pa-234m, Th-234 | U-238 | Chemically separated |
| Ba-137m | Cs-137 | (Ba-137m is always from Cs-137 in practice) |

Assign peaks to the long-lived parent (e.g., Ra-226 not Pb-214) unless there is clear evidence of chemical separation or disequilibrium.

## Source Validation

Before finalizing a source identification:

1. **Energy fit**: Source gamma energy must be within ~1 FWHM of the peak mean.
2. **Context check**: Use `source_info` on a candidate to check for common mis-IDs, associated sources, and whether it makes sense in context.
3. **Corroborating peaks**: If a source has multiple gammas, check if other peaks are detected. Consider:
   - Attenuation (especially at lower energies)
   - Low branching ratio lines may not produce detectable peaks
   - Nearby stronger peaks may obscure weaker lines
4. **Relative efficiency** (>2 peaks): Use `peak_based_relative_efficiency` on one nuclide at a time.
   - `residual_sigma` < 5 is good.
   - Higher values acceptable if: small peak next to large one, 511 keV annihilation interference, low-energy shielding effects, or secular equilibrium issues.
   - Try different `eqn_form` values (LnX, LnY, LnXLnY, FramEmpirical) if results are poor. Do not use `FramPhysicalModel` for validation.

## Escape and Sum Peaks

### Escape Peaks
When a gamma ray >1022 keV undergoes pair production in the detector:
- **Single escape peak**: At (gamma energy - 511) keV
- **Double escape peak**: At (gamma energy - 1022) keV

Use `escape_peak_check` to identify these. More common in HPGe detectors for high-energy gammas.

### Sum Peaks
Two types:
- **Cascade-sum**: Two gammas emitted simultaneously from same nuclide decay, detected together. Requires source close to detector (<~15 cm). Both contributing peaks have same nuclide source.
- **Random-sum (pile-up)**: Unrelated photons arriving simultaneously. Occurs at high count rates (>20% dead time for HPGe).

Use `sum_peak_check` to identify these.

### Pure Beta Emitters (Bremsstrahlung)
Pure beta emitters (P-32, Sr-90/Y-90, Tl-204, Sr-89, ...) produce no photopeaks - only a smooth
bremsstrahlung continuum that rises to a hump at ~60-150 keV and falls quasi-exponentially.
Use `beta_continuum_check` (requires a background spectrum loaded) - it background-subtracts,
excludes peak regions, and grades the net continuum shape.  The reported termination energy is a
statistics-limited LOWER bound on the beta endpoint.  X-ray-region peaks below ~120 keV do not
rule a beta source out (e.g. Tl-204's Hg K x-rays from its small EC branch).

## Detector Response Functions (DRF)

A DRF describes the detector's intrinsic efficiency as a function of energy. Required for:
- Activity fitting (`activity_fit`, `activity_fit_one_off`)
- MDA calculations (`currie_mda_calc` with nuclide/distance)
- Photopeak detection efficiency (`photopeak_detection_efficiency`)

Check if a DRF is loaded: `detector_efficiency_function_info` (no arguments).
List available DRFs: `available_detector_efficiency_functions`.
Load one: `load_detector_efficiency_function` by name, path, formula, or energy/efficiency pairs.

## Shielding and Geometry

### Geometry Types for Activity Fitting
- **Spherical**: Point source with spherical shielding. Default assumption.
- **CylinderSideOn**: Cylindrical source viewed from the side. For pipes, rods.
- **CylinderEndOn**: Cylindrical source viewed from the end. Also used for ground/surface contamination modeling.
- **Rectangular**: Rectangular source geometry.

### Source Types
- **Point**: Discrete source at center of shielding.
- **Trace**: Volumetric contamination distributed within a shielding layer. Specify via `trace_sources` in shielding definition.
- **Intrinsic**: Self-attenuating source (nuclide is part of the shielding material itself).

### Surface/Ground Contamination
Model using CylinderEndOn geometry with thin (~5 mm) "Air" shielding containing a trace source. Use large radius (~30 m) unless dimensions are specified.
