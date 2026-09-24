# Direct EPICS2023 photon migration validation

Historical baseline: `ab03e48`. Raw historical columns are read transiently and are not reproduced.

Absorption-edge neighborhoods include the full historical interpolation bracket around direct or historical discontinuities and are reported separately.

| Process | Away p95 | Away worst | Gate p95/worst | Edge p95/worst | Result |
|---|---:|---:|---:|---:|---|
| photoelectric | 1.1061% | 3.57556% | 1.25% / 4% | 0.104821% / 1.10319% | PASS |
| compton | 0.109713% | 0.172469% | 0.25% / 0.3% | 0.158971% / 0.191694% | PASS |
| rayleigh | 0.63573% | 6.03547% | 0.75% / 6.5% | 0.16198% / 0.533916% | PASS |
| pair_production | 0.0400808% | 0.0568256% | 0.05% / 0.1% | 0% / 0% | PASS |

## Contribution-weighted total attenuation

Away from edges: p95 0.712339%, worst 1.65817% (gate 0.5% / 2%): FAIL.

Edge neighborhoods (reported separately): p95 0%, worst 0%.
