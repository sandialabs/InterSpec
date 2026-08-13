# Bremsstrahlung migration review bundle

Open `review.html` in a browser. Select each element/material and incident
energy to compare the historical main/feature spectra, the migrated runtime
spectrum, and the deterministic Geant4 reference observable.

Review these merge-blocking questions:

- Is there any structured low-κ shift that is not understood?
- Is there an endpoint discontinuity or interpolation artifact?
- Is any material-dependent absolute-rate change unexplained?

The exhaustive Z=1–92 by 27-native-energy comparison against the Geant4
reference passes its numerical gates:

| Metric | Median | p95 | Worst |
|---|---:|---:|---:|
| Row RMS | 0.011927% | 0.054334% | 0.252031% |
| Integrated rate | 0.009052% | 0.049123% | 0.188068% |
| Mean κ | 0.009053% | 0.049147% | 0.188423% |

`summary.csv` contains absolute yields, rates, moments, and tail probabilities.
`spectra.csv` contains only 64-bin aggregate spectra/CDFs, not native G4
cross-section columns. `worst_20.csv` contains the worst exhaustive cases.
`metadata.json` records versions and fixture provenance.

After inspection, edit `manual_review.json`: set the three resolution fields to
`true`, record reviewer/time/notes, and set `status` to `approved`. The gate is:

```bash
python3 tools/prepare_cross_sections/check_provenance.py --require-manual-approval
```
