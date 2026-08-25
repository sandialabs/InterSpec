# GEANT4 Reference Data

This directory contains pre-computed GEANT4 reference efficiency values for the
standard benchmark configurations.  They let CeeLo be validated against GEANT4
without requiring a GEANT4 installation: `profiling/compare_validation.py` gates
these references against CeeLo's own results in `../ceelo_reference/`, and
`tests/test_cascade_summing.cpp` reads `cascade_summing_multi.csv` directly.

## How to regenerate reference data

1.  Build the GEANT4 validation harness:
    ```bash
    cmake -B build -S /path/to/mc_det_eff_plan \
          -DCEELO_BUILD_GEANT4_VALIDATION=ON \
          -DGeant4_DIR=/path/to/geant4/lib/cmake/Geant4
    cmake --build build -j8
    ```

2.  Export the geometry and macro from any example program, then run:
    ```bash
    ./build/tools/geant4_validation/ceelo_g4val \
        detector.gdml run.mac output.csv 4
    ```

3.  Copy `output.csv` to this directory with a descriptive name, e.g.:
    ```bash
    cp output.csv tests/data/geant4_reference/nai_3x3_10cm_multi.csv
    ```

4.  Commit the new reference file.

## File naming convention

`<detector>_<attenuators>_<source>_multi.csv` — one file per configuration,
holding every energy for that configuration (one row each).  For example
`nai_3x3_al1mm_pb2mm_15cm_multi.csv` is a 3"×3" NaI behind 1 mm Al + 2 mm Pb
with a point source at 15 cm.  Per-energy CSVs written by the harness are
assembled into one `_multi.csv` before committing.

## CSV format

```
# ceelo_g4val reference data
# geometry: detector.gdml
# macro: run.mac
energy_keV,fep_efficiency,fep_uncertainty,total_efficiency,total_uncertainty,num_events
662.0,0.01234,0.00012,0.03456,0.00020,1000000
```

## Current reference files

Generated with GEANT4 11.4.0, `FTFP_BERT` + `G4EmStandardPhysics_option4`.
Point-source configs are cone-biased at 64M events/energy; extended and shielded
sources are isotropic at 16M–32M, except configs 11 {2000, 3000} and 12 {662,
1173, 2614}, which were regenerated at 64M.  Per-file provenance notes are in the
CSV headers.

| File | Detector | Source | Config |
|------|----------|--------|--------|
| `nai_3x3_10cm_multi.csv` | 3"×3" NaI, bare | point, on-axis 10 cm | 1 |
| `nai_3x3_al1mm_10cm_multi.csv` | 3"×3" NaI, 1 mm Al endcap | point, on-axis 10 cm | 2 |
| `labr3_2x2_al05mm_5cm_multi.csv` | 2"×2" LaBr₃, 0.5 mm Al endcap | point, on-axis 5 cm | 3 |
| `czt_1x1x05cm_5cm_multi.csv` | 1×1×0.5 cm CZT box, bare | point, on-axis 5 cm | 5 |
| `nai_3x3_offaxis45_15cm_multi.csv` | 3"×3" NaI, bare | point, off-axis 45°, 15 cm | 6 |
| `nai_3x3_al1mm_pb2mm_15cm_multi.csv` | 3"×3" NaI, 1 mm Al + 2 mm Pb | point, on-axis 15 cm | 7 |
| `nai_3x3_al05mm_marinelli_water_multi.csv` | 3"×3" NaI, 0.5 mm Al | Marinelli beaker, water fill | 8 |
| `nai_3x3_fe05cm_shield_10cm_multi.csv` | 3"×3" NaI, bare | point, 10 cm, 0.5 cm Fe shield | 11 |
| `nai_3x3_ss304box_cellulose_15cm_multi.csv` | 3"×3" NaI, bare | SS304 box + cellulose, 15 cm | 12 |
| `hpge_gem35_coax_sharp_5cm_multi.csv` | GEM35-70 HPGe coax, sharp front edge | point, on-axis 5 cm | 25 |
| `hpge_gem35_coax_bullet_5cm_multi.csv` | GEM35-70 HPGe coax, bulletized edge + round-tipped bore | point, on-axis 5 cm | 26 |
| `cascade_summing_multi.csv` | 3"×3" NaI | Al cylinder, full-decay ion source | cascade-summing gate |

Configs 25 and 26 are a **matched pair**: the same crystal with a sharp and a
bulletized front edge. Each is gated against GEANT4 on its own, but the point of
the pair is the *difference* between them, which isolates the geometry change
from the physics residuals both configs share. They carry no dead layer, because
`write_gdml()` does not export one (see TODO.md), and their crystal is exported
as a single G4Polycone rather than boolean solids — a subtraction-based export
produced ~0.3% stuck-track warnings under cone bias, the polycone none.

The authoritative, up-to-date agreement figures are the tolerances encoded in
`profiling/compare_validation.py` and the validated-configuration table in
`DESIGN.md`; run

```bash
python3 profiling/compare_validation.py
```

to reproduce them.  The historical summary below is kept for context.

## MC vs G4 agreement summary (4M MC events vs 64M G4 events)

| Config | Max FEP discrepancy | Max total discrepancy | Notes |
|--------|--------------------|-----------------------|-------|
| 1 (NaI bare) | -0.96% | -0.40% | |
| 2 (NaI + Al) | -0.95% | -1.53% | |
| 3 (LaBr3 + Al) | -1.41% | -1.02% | |
| 5 (CZT) | -6.57% (1500 keV) | -0.77% | FEP < 1.5% for 30–800 keV |
| 6 (NaI off-axis) | -0.83% | -0.55% | |
| 7 (NaI + Al + Pb) | -0.99% (excl. 100 keV) | -1.10% | 100 keV has ~zero efficiency through Pb |
| 11 (NaI + Fe shield) | -0.6% | -1.2% | GDML uses rmin=1e-4 for sphere |
| 12 (NaI + SS304 box) | +2.1% | +1.9% | 16M G4, 7.6M MC |

## Known remaining physics discrepancies

| Effect | Energy range | Δ FEP | Root cause |
|--------|-------------|-------|------------|
| Doppler broadening | ≥1 MeV (all configs) | −0.5 to −1.0% | Missing bound-electron Compton profiles |
| Electron escape | ≥1 MeV (CZT only) | −3 to −7% | Thin detector; CSDA tracking vs G4 full e⁻ transport |

## Cascade true-coincidence summing reference

`cascade_summing_multi.csv` is the single source of truth for the cascade-summing gate. It is
read by **both** the C++ ctest (`tests/test_cascade_summing.cpp`) and
`profiling/compare_validation.py`, so the reference numbers are no longer hardcoded in the test.

Schema (one row per observable):
```
nuclide,estimator,name,lo,hi,peak_keV,g4,g4_sig,rlo,rhi,source_cm,emax
```
`g4,g4_sig` are the GEANT4 per-decay area over `[lo,hi]` keV and its 1σ; `rlo,rhi` is the
accepted engine/G4 ratio band; `source_cm` is the alcyl Al-cylinder radius/halfz (0.03 cm; 0.15
for Na-22). `estimator=full` gates a FullRealization summed-spectrum area; `estimator=cond`
gates the Conditional estimator (per-decay area = `A_FR·k_cond/k_full`) against the same G4 area.

Geometry ("alcyl"): 3"×3" NaI, source distributed in a solid Al cylinder (r = halfz = `source_cm`)
2 cm behind the front face. **The G4 GPS sampling volume MUST equal `source_cm`** — shrinking it
(e.g. a 0.024 cm stuck-track workaround) buries the decays deeper and adds L x-ray self-absorption
the engine does not model (≈+10% on total-L).

### How to regenerate

The G4 references come from `ceelo_g4val ... --histogram --correlated-gamma` at high stats
(8M decays Co-60/Ba-133/Co-57, 4M the rest) in the alcyl geometry; extract per-decay window
areas from the histogram and bake them into the rows. Example (Co-57):
```bash
# 1. export the alcyl GDML (the mm arg -> cm: pass 0.3 for the 0.03 cm cylinder; writes Co57.gdml)
./build/examples/cascade_g4_spectrum Co57 2000 2.0 0.3 260 alcyl
# 2. hand-author a GPS ion macro confining /gps/ion 27 57 in SrcMaterialPV (radius/halfz 0.03 cm).
#    CRITICAL: long-lived isotopes (Co-57 271 d, and Co-60/Ba-133/Am-241) do NOT decay without
#      /process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year
#    else "Events with deposit: 0". GPS volume MUST equal the source region (0.03), never shrunk.
source <your-geant4-install>/bin/geant4.sh
ceelo_g4val Co57.gdml run.mac co57_alcyl.csv --histogram --correlated-gamma --lowcut
# 3. per-decay area = window counts / (# Events); 1sigma = sqrt(counts) / (# Events).
```
(The efficiency columns in the harness CSV read 0 for an ion source — they key off a mono primary
energy — so the cascade gate uses the `--histogram` output, not those columns.)
The engine side (`our_cascade_multi.csv`) is produced by `build/examples/cascade_observables`.

## GEANT4 navigation notes

All cylindrical configs show ~0.06% "stuck track" warnings (`G4Exception: Stuck Track`) during
simulation. These arise from G4Tubs precision near r=0 when photons travel exactly on-axis.
GEANT4 resolves each automatically with a sub-micron push (`step = 1e-7 mm`); the effect on
efficiency is negligible and the rate is the same for bare-crystal and attenuated configs.
