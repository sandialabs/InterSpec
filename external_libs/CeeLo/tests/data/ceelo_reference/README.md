# CeeLo reference data read by tests

## `analytic_cascade_mc.csv`

The frozen Monte Carlo inputs of `tests/test_analytic_cascade_nuclides.cpp`
(`ANALYTIC_MC_REF`): the efficiency cache built by `compute()` and the FullRealization
reference summing factors for real nuclides. Freezing them keeps that test MC-free and
deterministic; it compares only the analytic summing combinatorics against the cached
FullRealization result.

To regenerate: set `REGENERATE_MC_VALUES` to 1 at the top of
`test_analytic_cascade_nuclides.cpp`, rebuild, run the test once (slow, real MC), set the
macro back to 0, and commit the updated CSV.

## CeeLo validation results (regenerated, not committed)

CeeLo's own full-energy-peak (FEP) and total efficiencies for the validated benchmark
configurations pair 1:1 with the GEANT4 references in
[`../geant4_reference/`](../geant4_reference/). They are regenerated from the current code
rather than committed, and `profiling/compare_validation.py` reads them from
`build/examples/` by default (`--mc-dir` to point elsewhere). It prints per-energy
FEP/total discrepancies and z-scores against the GEANT4 references.

`our_<config>_multi.csv` has one row per energy; columns
`energy_keV,fep_efficiency,fep_uncertainty,total_efficiency,total_uncertainty,num_events`
(the GEANT4 references' schema).

| File | Config | GEANT4 counterpart (`../geant4_reference/`) |
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
| `our_27_multi.csv` | 27 — GEM35-70 HPGe, bulletized, 2 cm + 0.5 cm Fe shell (0.75 keV window) | `hpge_gem35_bullet_fe05cm_2cm_multi.csv` |
| `our_28_multi.csv` | 28 — GEM35-70 HPGe, bulletized, 10 cm + 0.5 cm Fe shell (0.75 keV window) | `hpge_gem35_bullet_fe05cm_10cm_multi.csv` |
| `our_cascade_multi.csv` | cascade summing (6 nuclides, "alcyl" geom) | `cascade_summing_multi.csv` |

Configs 25 and 26 are a matched pair — the same crystal with a sharp and a bulletized
front edge — so the difference between their rows measures the bulletization effect
itself. See DESIGN.md → Validated Configurations for that comparison against GEANT4.

### Regenerating

Build with `-DCMAKE_BUILD_TYPE=Release` (or RelWithDebInfo); biasing is auto-enabled and
the precision target is ~0.3%. From `build/examples/`:

```bash
for c in 1 2 3 5 6 7 8 11 12 25 26; do
    ./benchmark_mc_configs --config $c --precision 0.003     # writes our_${c}_multi.csv (cwd)
done
for c in 27 28; do                                           # scored at their G4 refs' window
    ./benchmark_mc_configs --config $c --precision 0.003 --fep-window 0.75
done
./cascade_observables --out our_cascade_multi.csv            # 8M FullRealization / 4M Conditional
python3 ../../profiling/compare_validation.py
```

- Configs 27 and 28 **must** get `--fep-window 0.75`: their GEANT4 references were scored
  at a 0.75 keV half-window, every other config at the pinned 1.5 keV. Getting it wrong is
  silent: on config 27 the 1.5 keV window reads +2.6 / +3.5 ± 0.4% more FEP at 60 keV and
  +1.5 / +2.2 ± 0.4% at 88 keV, and 60/88 keV are documented skips, so the gate would not
  notice.
- `--precision p` targets p on **both** FEP and total, and the run continues until both are
  met; FEP is the binding one. Pass only `target_fep_rel_precision` for a cheaper run.
- Config 7 at 100 keV cannot reach 0.3% on FEP: it is `max_events`-limited at 200M events
  and ~4.5%. It is a documented skip in the gate.
- `our_cascade_multi.csv` holds per-decay summing observables
  (`nuclide,estimator,name,area,area_unc`) for the same observables the C++ ctest gate
  checks: `estimator=full` is a FullRealization summed-spectrum window area,
  `estimator=cond` the Conditional per-decay area (`A_FR·k_cond/k_full`).

Config 5 ≥1 MeV (CZT e⁻ escape) and config 7 at 100 keV (Pb K-edge) carry documented
tolerances/skips — see `profiling/compare_validation.py` and DESIGN.md "Validated
Configurations". Each CSV's `#` header records the run's event count and precision.
