# CeeLo Reference Results

CeeLo's own full-energy-peak (FEP) and total efficiency results for the validated
benchmark configurations, committed as a reproducible record **paired 1:1** with the
GEANT4 references in [`../geant4_reference/`](../geant4_reference/). Together they are the
"results from G4 and the current corresponding CeeLo result" that travel with the library.

`profiling/compare_validation.py` reads these by default and prints per-energy FEP/total
discrepancies and z-scores against the G4 references.

## Files

`our_<config>_multi.csv` — one row per energy; columns
`energy_keV,fep_efficiency,fep_uncertainty,total_efficiency,total_uncertainty,num_events`
(same schema as the G4 references).

| File | Config | G4 counterpart (`../geant4_reference/`) |
|------|--------|-----------------------------------------|
| `our_1_multi.csv`  | 1 — 3"×3" NaI, bare, 10 cm            | `nai_3x3_10cm_multi.csv` |
| `our_2_multi.csv`  | 2 — 3"×3" NaI, 1 mm Al, 10 cm         | `nai_3x3_al1mm_10cm_multi.csv` |
| `our_3_multi.csv`  | 3 — 2"×2" LaBr₃, 0.5 mm Al, 5 cm      | `labr3_2x2_al05mm_5cm_multi.csv` |
| `our_5_multi.csv`  | 5 — 1×1×0.5 cm CZT, bare, 5 cm        | `czt_1x1x05cm_5cm_multi.csv` |
| `our_6_multi.csv`  | 6 — 3"×3" NaI, off-axis 45°, 15 cm    | `nai_3x3_offaxis45_15cm_multi.csv` |
| `our_7_multi.csv`  | 7 — 3"×3" NaI, 1 mm Al + 2 mm Pb, 15 cm | `nai_3x3_al1mm_pb2mm_15cm_multi.csv` |
| `our_8_multi.csv`  | 8 — 3"×3" NaI, 0.5 mm Al, Marinelli (water) | `nai_3x3_al05mm_marinelli_water_multi.csv` |
| `our_11_multi.csv` | 11 — 3"×3" NaI, 0.5 cm Fe shield, 10 cm | `nai_3x3_fe05cm_shield_10cm_multi.csv` |
| `our_12_multi.csv` | 12 — 3"×3" NaI, SS304 box + cellulose | `nai_3x3_ss304box_cellulose_15cm_multi.csv` |
| `our_25_multi.csv` | 25 — GEM35-70 HPGe coax, sharp front edge, 5 cm | `hpge_gem35_coax_sharp_5cm_multi.csv` |
| `our_26_multi.csv` | 26 — GEM35-70 HPGe coax, bulletized edge + round-tipped bore, 5 cm | `hpge_gem35_coax_bullet_5cm_multi.csv` |
| `our_cascade_multi.csv` | cascade summing (6 nuclides, "alcyl" geom) | `cascade_summing_multi.csv` |

Configs 25 and 26 are a matched pair — the same crystal with a sharp and a
bulletized front edge — so the difference between their rows measures the
bulletization effect itself. See DESIGN.md → Validated Configurations for that
comparison against GEANT4.

### `our_cascade_multi.csv` (cascade true-coincidence summing)

Engine per-decay summing observables (`nuclide,estimator,name,area,area_unc`) for the same
observables the C++ ctest gate checks, paired with the G4 reference bands in
`../geant4_reference/cascade_summing_multi.csv`. `estimator=full` is a FullRealization
summed-spectrum window area; `estimator=cond` is the Conditional per-decay area
(`A_FR·k_cond/k_full`). Regenerate from `build/examples/` with:

```bash
./cascade_observables --out our_cascade_multi.csv   # 8M FullRealization / 4M Conditional
cp our_cascade_multi.csv ../../tests/data/ceelo_reference/
```

## How to regenerate / refresh

Built with `-DCMAKE_BUILD_TYPE=RelWithDebInfo` (or Release), auto-enabled biasing, precision
target ~0.3%. From `build/examples/`:

```bash
for c in 1 2 3 5 6 7 8 11 12 25 26 27 28; do
    ./benchmark_mc_configs --config $c --precision 0.003     # writes our_${c}_multi.csv (cwd)
done
cp our_*_multi.csv ../../tests/data/ceelo_reference/         # commit the refreshed snapshot
```

`--precision p` targets p on **both** FEP and total, and the run continues until **both** are
met. That was not always true: the two targets used to stop the run independently, so a run
ended at whichever converged first — always `total`, being the larger efficiency — and left FEP
short of the precision it was asked for. The set committed before Sep 7 2026 carries that:
FEP relative precision a median 1.44x worse than total, and worse than 2x on 29 of 122 rows.
Pass only the metric you care about (`target_fep_rel_precision` alone) if you want the cheaper
run; FEP is the binding one, so `total` comes along for free.

Config 7 at 100 keV cannot reach 0.3% on FEP: it is already `max_events`-limited at 200M events
and ~4.5%, and would need ~44 billion. It is a documented SKIP in the gate.

To compare a **fresh** dev run without overwriting this committed record:

```bash
python3 profiling/compare_validation.py --mc-dir build/examples
```

Config 5 ≥1 MeV (CZT e⁻ escape) and config 7 @100 keV (Pb K-edge) carry documented
tolerances/skips — see `profiling/compare_validation.py` and DESIGN.md "Validated
Configurations". Each CSV's `#` header records the run's event count / precision.
