# IDB planar HPGe uranium spectrum

`69_0044623U235.spe` is one unmodified spectrum, record 69, from the IAEA International Database of
Reference Gamma Spectra, version IDB-v2024-01, published by the IAEA Nuclear Data Section:

- Dataset: https://doi.org/10.61092/iaea.gaef-dhff
- Data descriptor: https://doi.org/10.1038/s41597-025-05455-1
- Licence: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
- Data provider (per the IDB metadata): United States Support Program

The file keeps its original IDB name and is unmodified.

## Measurement

| | |
|---|---|
| Material | NBS SRM 969 uranium enrichment standard, U3O8 powder |
| Certified (wt%, 2 sigma) | U-234 0.0359 ± 0.0068, U-235 4.4623 ± 0.0032, U-236 0.0068 ± 0.0002, U-238 95.495 ± 0.0032 |
| Certificate | Reference date 1986-01-01; https://nds.iaea.org/media-idb/certificates/LEU_SRM969.pdf#page=103 |
| Detector | Planar HPGe, FWHM 0.65 keV, 4096 channels at 0.075 keV/channel (0-307 keV) |
| Acquired | 2002-02-28 23:39, 1599 s live / 1800 s real |

Decay over the 16 years between the certificate and the measurement changes these fractions by much
less than their uncertainties, so the certified values are the truth.

## Use

`test_LlmIsotopicsTool.cpp` runs the shipped "HPGe U (63-205 keV)" preset on it, with the
"2.5x1.5 Planar HPGe" detector response.

This record was chosen from the 17 repeat measurements of this standard on this detector because that
preset's result for it is typical of them: in a development study, U-235 came out
0.08% below the certified value (i.e. 0.9992*certified_value), next to the median of the 17, and the fit's chi2/dof was typical too.
