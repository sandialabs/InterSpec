# EPICS2023 process-specific photon table report

Every native source node in 10 keV–20 MeV is validated, but only endpoints, thresholds, duplicated-edge flank/edge pairs, and adaptively required nodes are retained. The reported errors include uint16 value decoding and 16 logarithmic probes per interval.

The form-factor curve is generator-only: its error columns validate the curve used to build the offline Rayleigh CDF, while its runtime retained-node and byte counts are zero.

The shared energy-coordinate pool stores 5,418 float32 values (21,672 bytes); process grid columns below are uint16 pool indices.

Rayleigh total cross sections use 25 shared adjacent-Z adaptive grids containing 9,566 total float32 coordinates (38,264 bytes).

| Process | Native | Retained | Coeff bytes | Grid bytes | Nodes min/med/p95/max | p95/max error | Max integrated |
|---|---:|---:|---:|---:|---:|---:|---:|
| rayleigh | 10,413 | 37,604 | 75,208 | 38,264 | 317/386/415/418 | 0.06274%/0.09996% | 0.02605% |
| compton | 8,288 | 3,693 | 7,386 | 7,386 | 35/37/40/43 | 0.07406%/0.09982% | 0.01902% |
| pair_production | 38,505 | 21,983 | 43,966 | 43,966 | 192/225/239/243 | 0.03714%/0.05000% | 0.01963% |
| photoelectric | 406,218 | 5,588 | 11,176 | 11,176 | 43/57/65/70 | 0.07682%/0.09999% | 0.04872% |
| k_photoelectric | 25,061 | 4,677 | 9,354 | 9,354 | 38/48/55/60 | 0.08337%/0.09999% | 0.04925% |
| form_factor | 106,810 | 0 | 0 | 0 | 0/0/0/0 | 0.07359%/0.10000% | 0.03188% |
| scattering_function | 17,376 | 9,201 | 18,402 | 36,804 | 75/95/102/103 | 0.08148%/0.10000% | 0.00602% |

Angular distribution gates:
maximum CDF error = 0.038119% (limit 0.25%); maximum mean-cosine error = 0.0002151 (limit 0.002).
Offline Rayleigh inverse-CDF uint16 quantization error = 7.63e-06 (limit 8e-06).
Runtime Rayleigh 256-node inverse-CDF discretization error = 0.128618% (limit 0.25%).
