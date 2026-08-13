# ESTAR stopping migration validation

Historical baseline: `ab03e48`.

| Quantity | Median | p95 | Worst | Worst case | Gate | Result |
|---|---:|---:|---:|---|---|---|
| collision_stopping | 0.164061% | 2.96347% | 15.8202% | Z=18, energy_keV=20000.0 | p95 ≤1%; worst ≤3% | FAIL |
| collision_stopping_historical_anchor_elements | 0.0108232% | 0.0777324% | 0.128635% | Z=29, energy_keV=150.0 | diagnostic | DIAGNOSTIC |
| csda_range | 0.337926% | 3.44405% | 10.7532% | Z=92, energy_keV=10.0 | worst ≤1% | FAIL |
| deterministic_deposited_energy | 0.211933% | 2.8573% | 15.3713% | Z=35, energy_keV=20000.0, old_range_fraction=0.1 | worst ≤1% | FAIL |

Overall: FAIL.

Direct ESTAR agrees closely at the 12 historical correction-surface anchors. The failed historical-result gate is concentrated in elements that the old runtime approximated by interpolation between those anchors, especially across elemental gas/solid density-effect changes. The direct 92×53 source table and its float32-log representation pass their separate source-fidelity checks.
