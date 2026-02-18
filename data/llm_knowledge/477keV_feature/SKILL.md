---
name: 477keV-feature
description: Background on common reaction channels that can contribute to or be discussed with the ~477 keV spectral feature, especially lithium and boron reaction pathways.
metadata:
  maturity: draft-example
  verification: Not robustly developed or independently checked.
---

# 477 keV Feature: Lithium and Boron Reaction Context

This skill provides qualitative context for discussion of a feature near 477 keV, with emphasis on lithium and boron reactions that may be relevant in neutron fields or neutron-interaction environments.

## Relevant channels to consider

- `Li(a,a)`: Elastic alpha scattering on lithium produces a strong Doppler Broadened (triangular) feature, that especially for HPGe detectors, may not fit a peak.
- `Li(n,n)`: Elastic/inelastic neutron scattering on lithium can create at ~478 keV.
- `B(n,n)`: Neutron scattering on boron can create at ~478 keV.
- `B(n,g)`: Neutron capture on boron (specifically B-10) can produce prompt gamma emission, resulting in a doppler broadened peak (10 to 15 keV FWHM) at ~478 keV.
- `B(n,a)`: Can also produce a signal at ~478 keV.
- `Be7`: produces a lone gamma at 477.6 keV.
- Other multi-gamma nuclides produce gammas near this region

## Feature detection

The Doppler Broadening (triangle) shape may mean that a peak will not be fit for that feature.  A work-around to detect the elevated region is to make three tool calls to get the counts in equally large energy ranges below, around, and above the range in question, i.e.,:
- `get_counts_in_energy_range( {"lowerEnergy": 458, "upperEnergy": 468} )`
- `get_counts_in_energy_range( {"lowerEnergy": 473, "upperEnergy": 483} )`
- `get_counts_in_energy_range( {"lowerEnergy": 488, "upperEnergy": 498} )`
If the center range is substantially elevated over both side regions, there is likely a feature there.

