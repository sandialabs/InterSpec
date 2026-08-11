# Multi-curve test spectra

Deterministic simulated two-stacked-uranium-disk spectra and backgrounds from the 2026-07
RelActCalcAuto multi-curve review harness (`scratch/multicurve_2026-07/`, `MulticurveStudy simulate`,
seed 20260725).  Geometry: two U disks r = 2 cm (front thickness varies, back 3 mm), both behind a
1.0 g/cm2 Fe case, DRF "ORTEC Detective-X_LANL_100cm" (data/common_drfs.tsv), FWHM(E) =
0.566 + 0.0393*sqrt(E) keV, pure-Gaussian peaks, GADRAS Detective-X background continuum.

| file | scenario |
|---|---|
| `easy_two_disk_fore.n42` | "easy": front disk 1 mm DU (0.2 % U235), back 3 mm HEU (93.3 %), 4 h-equivalent at 1 m (live time 14400 s, act_scale 5) |
| `easy_background_long.n42` | matching 3x-long background |

Truth values (enrichments, areal densities, per-line expected counts) are in the harness
`runs/*/truth.txt`; the values asserted in `test_RelActCalcAuto_MultiCurve.cpp` are hardcoded from there.
