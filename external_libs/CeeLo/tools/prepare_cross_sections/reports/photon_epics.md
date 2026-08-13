# EPICS2023 process-specific photon table report

Every native source node in 10 keV–10 MeV is validated, but only endpoints, thresholds, duplicated-edge flank/edge pairs, and adaptively required nodes are retained. The reported errors include uint16 value decoding and 16 logarithmic probes per interval.

The form-factor curve is generator-only: its error columns validate the curve used to build the offline Rayleigh CDF, while its runtime retained-node and byte counts are zero.

The shared energy-coordinate pool stores 5,040 float32 values (20,160 bytes); process grid columns below are uint16 pool indices.

Rayleigh total cross sections use 23 shared adjacent-Z adaptive grids containing 8,156 total float32 coordinates (32,624 bytes).

| Process | Native | Retained | Coeff bytes | Grid bytes | Nodes min/med/p95/max | p95/max error | Max integrated |
|---|---:|---:|---:|---:|---:|---:|---:|
| rayleigh | 9,470 | 32,624 | 65,248 | 32,624 | 283/357/388/392 | 0.05917%/0.09991% | 0.02456% |
| compton | 6,674 | 3,125 | 6,250 | 6,250 | 32/34/35/36 | 0.07352%/0.09998% | 0.02283% |
| pair_production | 35,494 | 20,518 | 41,036 | 41,036 | 184/224/238/241 | 0.03701%/0.05000% | 0.01745% |
| photoelectric | 352,432 | 4,877 | 9,754 | 9,754 | 42/53/60/68 | 0.07792%/0.09996% | 0.04737% |
| k_photoelectric | 22,501 | 4,145 | 8,290 | 8,290 | 35/46/52/58 | 0.08432%/0.09999% | 0.04923% |
| form_factor | 101,209 | 0 | 0 | 0 | 0/0/0/0 | 0.07346%/0.10000% | 0.03188% |
| scattering_function | 16,382 | 8,602 | 17,204 | 34,408 | 75/94/102/103 | 0.08180%/0.10000% | 0.00602% |

Angular distribution gates:
maximum CDF error = 0.038119% (limit 0.25%); maximum mean-cosine error = 0.0002151 (limit 0.002).
Offline Rayleigh inverse-CDF uint16 quantization error = 7.63e-06 (limit 8e-06).
Runtime Rayleigh 256-node inverse-CDF discretization error = 0.128618% (limit 0.25%).
