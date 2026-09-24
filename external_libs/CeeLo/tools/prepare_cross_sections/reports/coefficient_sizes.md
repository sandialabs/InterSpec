# Cross-section coefficient and resident-data accounting

Shared grids are listed separately. Size increases are accepted only where the committed accuracy report passes its gate.

| Baseline | Process | Before coeffs / encoding | After retained / native | Before bytes | After bytes | Change | p95 / max error |
|---|---|---:|---:|---:|---:|---:|---:|
| main | rayleigh | 4,578 / float32 | 37,604 / 10,413 | 18,312 | 114,360 | +96,048 (+524.5%) | 0.062743% / 0.0999626% |
| main | compton | 4,578 / float32 | 3,693 / 8,288 | 18,312 | 15,948 | -2,364 (-12.9%) | 0.0740637% / 0.0998173% |
| main | pair_production | 4,578 / float32 | 21,983 / 38,505 | 18,312 | 89,108 | +70,796 (+386.6%) | 0.037135% / 0.0499998% |
| main | photoelectric | 4,578 / float32 | 5,588 / 406,218 | 18,312 | 23,528 | +5,216 (+28.5%) | 0.0768154% / 0.0999878% |
| main | k_photoelectric | 4,578 / float32 | 4,677 / 25,061 | 18,312 | 19,884 | +1,572 (+8.6%) | 0.083375% / 0.0999879% |
| main | form_factor | 1,656 / float32 | 0 / 106,810 | 6,624 | 0 | -6,624 (-100.0%) | 0.0735897% / 0.0999979% |
| main | scattering_function | 1,656 / float32 | 9,201 / 17,376 | 6,624 | 56,382 | +49,758 (+751.2%) | 0.0814771% / 0.0999962% |
| main | bremsstrahlung | 79,488 / float32 | 79,488 / 79,488 | 318,192 | 159,584 | -158,608 (-49.8%) | 0.010051% / 0.243092% |
| main | collision_stopping | 180 / float32 correction ratio | 4,876 / 4,876 | 720 | 10,488 | +9,768 (+1356.7%) | 0.00165736% / 0.00194615% |
| main | radiative_stopping | 180 / float32 correction ratio | 4,876 / 4,876 | 720 | 10,488 | +9,768 (+1356.7%) | 0.00340054% / 0.00400204% |
| main | K_relaxation | 1,288 / float32 | 2,156 / 2,154 | 5,520 | 11,024 | +5,504 (+99.7%) | 3.397e-06% / 5.30243e-06% |
| main | L_relaxation | 4,300 / float32 | 3,090 / 3,083 | 18,400 | 19,560 | +1,160 (+6.3%) | 3.397e-06% / 5.30243e-06% |
| feature/FullSpecApprox | rayleigh | 4,578 / float32 | 37,604 / 10,413 | 18,312 | 114,360 | +96,048 (+524.5%) | 0.062743% / 0.0999626% |
| feature/FullSpecApprox | compton | 4,578 / float32 | 3,693 / 8,288 | 18,312 | 15,948 | -2,364 (-12.9%) | 0.0740637% / 0.0998173% |
| feature/FullSpecApprox | pair_production | 4,578 / float32 | 21,983 / 38,505 | 18,312 | 89,108 | +70,796 (+386.6%) | 0.037135% / 0.0499998% |
| feature/FullSpecApprox | photoelectric | 4,578 / float32 | 5,588 / 406,218 | 18,312 | 23,528 | +5,216 (+28.5%) | 0.0768154% / 0.0999878% |
| feature/FullSpecApprox | k_photoelectric | 4,578 / float32 | 4,677 / 25,061 | 18,312 | 19,884 | +1,572 (+8.6%) | 0.083375% / 0.0999879% |
| feature/FullSpecApprox | form_factor | 1,656 / float32 | 0 / 106,810 | 6,624 | 0 | -6,624 (-100.0%) | 0.0735897% / 0.0999979% |
| feature/FullSpecApprox | scattering_function | 1,656 / float32 | 9,201 / 17,376 | 6,624 | 56,382 | +49,758 (+751.2%) | 0.0814771% / 0.0999962% |
| feature/FullSpecApprox | bremsstrahlung | 76,544 / uint16 | 79,488 / 79,488 | 153,692 | 159,584 | +5,892 (+3.8%) | 0.010051% / 0.243092% |
| feature/FullSpecApprox | collision_stopping | 180 / float32 correction ratio | 4,876 / 4,876 | 720 | 10,488 | +9,768 (+1356.7%) | 0.00165736% / 0.00194615% |
| feature/FullSpecApprox | radiative_stopping | 180 / float32 correction ratio | 4,876 / 4,876 | 720 | 10,488 | +9,768 (+1356.7%) | 0.00340054% / 0.00400204% |
| feature/FullSpecApprox | K_relaxation | 1,288 / float32 | 2,156 / 2,154 | 5,520 | 11,024 | +5,504 (+99.7%) | 3.397e-06% / 5.30243e-06% |
| feature/FullSpecApprox | L_relaxation | 4,300 / float32 | 3,090 / 3,083 | 18,400 | 19,560 | +1,160 (+6.3%) | 3.397e-06% / 5.30243e-06% |
| feature/env-scatter | rayleigh | 4,578 / float32 | 37,604 / 10,413 | 18,312 | 114,360 | +96,048 (+524.5%) | 0.062743% / 0.0999626% |
| feature/env-scatter | compton | 4,578 / float32 | 3,693 / 8,288 | 18,312 | 15,948 | -2,364 (-12.9%) | 0.0740637% / 0.0998173% |
| feature/env-scatter | pair_production | 4,578 / float32 | 21,983 / 38,505 | 18,312 | 89,108 | +70,796 (+386.6%) | 0.037135% / 0.0499998% |
| feature/env-scatter | photoelectric | 4,578 / float32 | 5,588 / 406,218 | 18,312 | 23,528 | +5,216 (+28.5%) | 0.0768154% / 0.0999878% |
| feature/env-scatter | k_photoelectric | 4,578 / float32 | 4,677 / 25,061 | 18,312 | 19,884 | +1,572 (+8.6%) | 0.083375% / 0.0999879% |
| feature/env-scatter | form_factor | 1,656 / float32 | 0 / 106,810 | 6,624 | 0 | -6,624 (-100.0%) | 0.0735897% / 0.0999979% |
| feature/env-scatter | scattering_function | 1,656 / float32 | 9,201 / 17,376 | 6,624 | 56,382 | +49,758 (+751.2%) | 0.0814771% / 0.0999962% |
| feature/env-scatter | bremsstrahlung | 76,544 / uint16 | 79,488 / 79,488 | 153,692 | 159,584 | +5,892 (+3.8%) | 0.010051% / 0.243092% |
| feature/env-scatter | collision_stopping | 180 / float32 correction ratio | 4,876 / 4,876 | 720 | 10,488 | +9,768 (+1356.7%) | 0.00165736% / 0.00194615% |
| feature/env-scatter | radiative_stopping | 180 / float32 correction ratio | 4,876 / 4,876 | 720 | 10,488 | +9,768 (+1356.7%) | 0.00340054% / 0.00400204% |
| feature/env-scatter | K_relaxation | 1,288 / float32 | 2,156 / 2,154 | 5,520 | 11,024 | +5,504 (+99.7%) | 3.397e-06% / 5.30243e-06% |
| feature/env-scatter | L_relaxation | 4,300 / float32 | 3,090 / 3,083 | 18,400 | 19,560 | +1,160 (+6.3%) | 3.397e-06% / 5.30243e-06% |

## Shared grids

| Baseline | Grid | Before bytes | After bytes | Note |
|---|---|---:|---:|---|
| main | historical_photon_energy_grids | 18,312 | 0 | Previously shared by five photon processes; never double-counted above. |
| main | historical_scattering_function_grid | 72 | 0 | Direct EPDL now uses per-element process-specific grids. |
| main | historical_form_factor_grid | 72 | 0 | Direct EPDL now uses per-element process-specific grids. |
| main | ESTAR_energy_grid | 0 | 212 | Shared by collision and radiative stopping; counted once. |
| main | EPDL_photon_energy_coordinate_pool | 0 | 21,672 | Bit-exact float32 pool shared by five process-specific uint16 index streams. |
| main | Rayleigh_offline_sampling_CDF | 0 | 51,200 | Normalized uint16 CDF for all elements plus one shared float32 x grid; generated offline. |
| feature/FullSpecApprox | historical_photon_energy_grids | 18,312 | 0 | Previously shared by five photon processes; never double-counted above. |
| feature/FullSpecApprox | historical_scattering_function_grid | 72 | 0 | Direct EPDL now uses per-element process-specific grids. |
| feature/FullSpecApprox | historical_form_factor_grid | 72 | 0 | Direct EPDL now uses per-element process-specific grids. |
| feature/FullSpecApprox | ESTAR_energy_grid | 0 | 212 | Shared by collision and radiative stopping; counted once. |
| feature/FullSpecApprox | EPDL_photon_energy_coordinate_pool | 0 | 21,672 | Bit-exact float32 pool shared by five process-specific uint16 index streams. |
| feature/FullSpecApprox | Rayleigh_offline_sampling_CDF | 0 | 51,200 | Normalized uint16 CDF for all elements plus one shared float32 x grid; generated offline. |
| feature/env-scatter | historical_photon_energy_grids | 18,312 | 0 | Previously shared by five photon processes; never double-counted above. |
| feature/env-scatter | historical_scattering_function_grid | 72 | 0 | Direct EPDL now uses per-element process-specific grids. |
| feature/env-scatter | historical_form_factor_grid | 72 | 0 | Direct EPDL now uses per-element process-specific grids. |
| feature/env-scatter | ESTAR_energy_grid | 0 | 212 | Shared by collision and radiative stopping; counted once. |
| feature/env-scatter | EPDL_photon_energy_coordinate_pool | 0 | 21,672 | Bit-exact float32 pool shared by five process-specific uint16 index streams. |
| feature/env-scatter | Rayleigh_offline_sampling_CDF | 0 | 51,200 | Normalized uint16 CDF for all elements plus one shared float32 x grid; generated offline. |

## Measured compiled resident-data sections

| Object | Bytes |
|---|---:|
| `photon_epics_data.cpp.o` | 392,084 |
| `element_data.cpp.o` | 183,528 |
| `relaxation_epics_data.cpp.o` | 29,384 |
| `estar_stopping_data.cpp.o` | 21,188 |

Compiled resident data total: 626,184 bytes.
Initialized singleton + one thread TLS: 797 bytes.
Combined for one active thread: 626,981 bytes (0.598 MiB).
