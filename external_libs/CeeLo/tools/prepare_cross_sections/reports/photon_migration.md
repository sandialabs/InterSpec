# Direct EPICS2023 photon migration validation

Historical baseline: `ab03e48`. Raw historical columns are read transiently and are not reproduced.

Absorption-edge neighborhoods include the full historical interpolation bracket around direct or historical discontinuities and are reported separately.

| Process | Away p95 | Away worst | Gate p95/worst | Edge p95/worst | Result |
|---|---:|---:|---:|---:|---|
| photoelectric | 1.10645% | 3.58305% | 1.25% / 4% | 0.106351% / 1.07092% | PASS |
| compton | 0.109545% | 0.158785% | 0.25% / 0.3% | 0.16054% / 0.192947% | PASS |
| rayleigh | 0.634881% | 6.03648% | 0.75% / 6.5% | 0.165858% / 0.522648% | PASS |
| pair_production | 0.0412566% | 0.0637807% | 0.05% / 0.1% | 0% / 0% | PASS |

## Contribution-weighted total attenuation

Away from edges: p95 0.717277%, worst 1.64468% (gate 0.5% / 2%): FAIL.

Edge neighborhoods (reported separately): p95 0%, worst 0%.
