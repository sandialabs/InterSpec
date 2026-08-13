# Third-party notices and data attribution

CeeLo itself is distributed under the GNU Lesser General Public License, version
2.1 — see [`LICENSE.txt`](LICENSE.txt).

This file records everything CeeLo depends on, builds against, or derives data
from. It is split into three parts:

1. **Software dependencies** — libraries CeeLo links or builds against.
   None are vendored into this repository.
2. **Evaluated data** — the numerical tables compiled into the library, and
   where each one comes from.
3. **Published methods** — algorithms implemented from the literature.

Nothing in `src/` is copied from another code base. Where an implementation
follows a published method that other codes also implement, the file says so at
the point of use; see in particular
[`src/transport/ComptonScatter.cpp`](src/transport/ComptonScatter.cpp),
[`src/io/Pchip.h`](src/io/Pchip.h), and
[`src/cascade/BetaSpectrum.cpp`](src/cascade/BetaSpectrum.cpp).

---

## 1. Software dependencies

None of these are redistributed with CeeLo; all are resolved at build time.

| Dependency | Licence | Required? | Used by |
|---|---|---|---|
| [Eigen](https://eigen.tuxfamily.org) 3 | MPL 2.0 | **Yes** | Vector/matrix algebra throughout `src/` |
| [Boost.Test](https://www.boost.org) | Boost Software License 1.0 | Tests only | `tests/` |
| [SandiaDecay](https://github.com/sandialabs/SandiaDecay) | LGPL 2.1 (NTESS) | Optional (`CEELO_WITH_SANDIADECAY`, default OFF) | `src/cascade/SandiaDecayCascade.*` — nuclear decay data for cascade summing |
| [RapidXml](http://rapidxml.sourceforge.net) | Boost Software License 1.0 / MIT | Optional (`CEELO_BUILD_RESPONSE`, default OFF) | `src/io/DetectorResponse.cpp` XML payload. Deliberately **not** vendored: the host application's copy is used, for ODR safety |
| [GEANT4](https://geant4.web.cern.ch) 11.4.0 | [Geant4 Software License](https://geant4.web.cern.ch/download/license) | Optional, developer-only | `tools/geant4_validation/` produces simulation and aggregate-observable validation fixtures. No production generator extracts GEANT4 coefficient tables. The CeeLo **library** has no GEANT4 dependency |
| [xraylib](https://github.com/tschoonj/xraylib) 4.2.1 | BSD 3-clause | Developer-only (data regeneration) | `tools/prepare_cross_sections/generate_element_support.py`; retained licence: `tools/prepare_cross_sections/source_notices/xraylib_license_all.txt` |

GEANT4 is additionally the reference implementation CeeLo is validated against.
Please cite it as:

- S. Agostinelli et al., "Geant4 — a simulation toolkit",
  *Nucl. Instrum. Meth. A* **506** (2003) 250–303.
- J. Allison et al., "Geant4 developments and applications",
  *IEEE Trans. Nucl. Sci.* **53** (2006) 270–278.
- J. Allison et al., "Recent developments in Geant4",
  *Nucl. Instrum. Meth. A* **835** (2016) 186–225.

---

## 2. Evaluated data compiled into the library

CeeLo ships all of its physics data as compiled-in C++ tables; there are no
runtime data files. Each table below carries a provenance banner in its own
source file.

The machine-readable source lock at
[`tools/prepare_cross_sections/sources.lock.json`](tools/prepare_cross_sections/sources.lock.json)
records exact URLs, versions, hashes, access dates, terms, citations, and
transformations for every generated table. Upstream archives are preparation
inputs only and are not committed.

### 2.1 Photon cross-sections — `src/cross_sections/photon_epics_data.cpp`

Photoelectric, incoherent (Compton), coherent (Rayleigh) and pair-conversion
cross-sections for Z = 1–92, on per-element adaptive energy grids.

- **Magnitudes and angular factors**: direct **EPICS2023 EPDL** records:
  MF=23/MT=502 Rayleigh, MT=504 Compton, MT=516 total pair production,
  MT=522 total photoelectric, MT=534 K-shell photoelectric, and MF=27/MT=502
  and 504 form factors and scattering functions.
- **Evaluation lineage**: EPICS / EPDL, the Evaluated Photon Data Library —
  D.E. Cullen, J.H. Hubbell and L. Kissel, *EPDL97: the Evaluated Photon Data
  Library, '97 Version*, UCRL-50400 Vol. 6 Rev. 5 (LLNL, 1997), and the
  IAEA Nuclear Data Services EPICS2023 release.
- Source nodes, thresholds, and both sides of discontinuities are retained on
  process-specific grids. Adaptive probes enforce the committed representation
  and integrated-error gates; see `reports/photon_epics.*`.
- **Compton subshell profiles** (occupancy, binding energy, per-electron J(0)):
  Biggs Compton profiles, via BSD-licensed xraylib 4.2.1/EPDL.

### 2.2 Bremsstrahlung — Seltzer-Berger tables in `element_data.cpp`

χ(Z, T, κ) on a 32 κ × 27 native electron-energy grid from 10 keV to 20 MeV,
Z=1–92. The source is the **NIST EPQ Berger–Seltzer table set** pinned at commit
`ce1e5895892b17854d857e596250e4ae70d6b274`. NIST-authored EPQ is a U.S.
Government work; its full public-domain notice and disclaimer are retained at
[`source_notices/NIST_EPQ_LicenseFile.txt`](tools/prepare_cross_sections/source_notices/NIST_EPQ_LicenseFile.txt).
The underlying reference is S.M. Seltzer and M.J. Berger, *Atomic Data and
Nuclear Data Tables* **35** (1986) 345–418.

The generator uniformly constructs κ=0.025 by quadratic interpolation through
EPQ κ={0, 0.05, 0.1} and κ=0.075 by cubic interpolation through
κ={0, 0.05, 0.1, 0.15}. It copies every other equal-κ value, omits EPQ κ=0.125,
does no per-element fitting, applies no row normalization, and emits `float32`.
EPQ `phi_rad` is used only as an independent first-moment checksum.

PENELOPE-2008 describes its related 32-node database as supplied by S. Seltzer,
but no independently redistributable copy of that later low-κ revision was
used. The fixed transformation above was validated against aggregate/numerical
Geant4 behavior; it is not a transcription, fitted reproduction, or coefficient
extraction from G4EMLOW. No G4EMLOW coefficient payload is present.

### 2.3 Atomic relaxation — `src/cross_sections/relaxation_epics_data.cpp`

K-shell fluorescence yields and line energies/intensities; per-subshell L1/L2/L3
yields, lines, and Coster-Kronig transfer yields. The relaxation-only table
extends through Z=99 to cover radioactive-decay daughters such as Np-237;
photon/electron transport remains limited to Z=92.

- Direct **EPICS2023 EADL** MF=28/MT=533, the Evaluated Atomic Data Library —
  S.T. Perkins et al., UCRL-50400 Vol. 30 Rev. 1 (LLNL, 1997), distributed by
  IAEA Nuclear Data Services. No Geant4 extraction layer is used.
- Cross-checked against W. Bambynek et al., *Rev. Mod. Phys.* **44** (1972) 716;
  M.O. Krause, *J. Phys. Chem. Ref. Data* **8** (1979) 307; and J.H. Hubbell
  et al., *J. Phys. Chem. Ref. Data* **23** (1994) 339.

### 2.4 Electron stopping powers — `src/physics/estar_stopping_data.cpp`

- Collision and radiative mass stopping powers: **NIST ESTAR** (M.J. Berger et
  al., NIST Standard Reference Database 124; M.J. Berger and S.M. Seltzer,
  NBSIR 82-2550), generated for every Z=1–92 at all 53 standard nodes from
  10 keV through 20 MeV. Runtime interpolation is linear in log energy/log
  stopping power. There is no representative-element correction surface.
- The downloadable NIST STAR program/database establishes algorithm and
  database provenance. NIST's download page permits redistribution and
  modification when its notices are retained and changes are identified; the
  source lock records the archive and component hashes.
- Mean excitation energies retained for lower-level models: **ICRU Report 49**
  (1993), Table 5.1.

### 2.5 Material compositions — `src/materials/BuiltinMaterials.cpp`

- Dry air at sea level: **ICRU Report 37** / NIST standard composition.
- Soil and other reference compositions: **IAEA-TECDOC-1011** (1998), Table 1,
  and the CRC Handbook of Chemistry and Physics.
- Atomic weights: IUPAC standard atomic weights (via xraylib; IUPAC 2021
  fallback table).

### 2.6 Nuclear decay data (optional cascade-summing path)

Cascade summing corrections use **SandiaDecay**, whose `sandia.decay.xml` is
derived from the **ENSDF** (Evaluated Nuclear Structure Data File) datasheets
maintained by NNDC/BNL, with manual corrections applied. SandiaDecay is not
vendored here and the cascade adapter is off by default. See
<https://www.nndc.bnl.gov/ensdf/>.

### 2.7 GEANT4 reference results — validation fixtures

Efficiency fixtures under `tests/data/geant4_reference/` and aggregate spectra,
rates, CDFs, and histograms under `validation/` are results produced by running
or evaluating configured GEANT4 models. They are project validation
measurements, not raw GEANT4 cross-section columns. Each fixture directory
records its version, configuration, and statistical provenance. Production
generators never consume these results as fitted correction tables.

---

## 3. Published methods implemented from the literature

Implemented from the cited publications, not ported from any implementation.

| Method | Where | Reference |
|---|---|---|
| Klein-Nishina angle sampling | `src/transport/ComptonScatter.cpp` | Butcher & Messel, *Nucl. Phys.* **20** (1960) 15 |
| Doppler broadening (impulse approximation, analytic subshell profiles) | `src/transport/ComptonScatter.cpp` | Ribberfors, *Phys. Rev. B* **12** (1975) 2067; PENELOPE's factorization |
| Photoelectron angular distribution | `src/physics/ElectronCsda.h` | Sauter (1931) relativistic dipole formula |
| Pair-production opening angle | `src/transport/PhotonTransport.cpp` | Koch & Motz, *Rev. Mod. Phys.* **31** (1959) 920 (modified Tsai form) |
| Multiple scattering (condensed history) | `src/physics/ElectronCsda.cpp` | Highland, *Nucl. Instrum. Meth.* **129** (1975) 497; Molière; Goudsmit & Saunderson, *Phys. Rev.* **57** (1940) 24. Background: Nelson et al., EGS4 (SLAC-265, 1985); Ivanchenko et al., GEANT4 EM physics (2010) |
| Radiation length | `src/physics/ElectronCsda.cpp` | Tsai, *Rev. Mod. Phys.* **46** (1974) 815 |
| Monotone cubic interpolation (PCHIP) | `src/io/Pchip.cpp` | Fritsch & Carlson, *SIAM J. Numer. Anal.* **17** (1980) 238. Derivative rule follows `scipy.interpolate.PchipInterpolator` (SciPy, BSD 3-clause) |
| log Γ(z) via Lanczos approximation | `src/cascade/BetaSpectrum.cpp` | Lanczos, *J. SIAM Numer. Anal. Ser. B* **1** (1964) 86. The g=7, n=9 coefficients are Godfrey's openly published set (also in Boost.Math) |
| 32-point Gauss-Legendre quadrature | `src/transport/ComptonScatter.cpp` | Abramowitz & Stegun, *Handbook of Mathematical Functions* (NBS, 1964), Table 25.4 — public domain |
| Efficiency transfer | `src/io/EfficiencyTransfer.cpp` | EFFTRAN: Vidmar, *Nucl. Instrum. Meth. A* **550** (2005) 603 |

---

## Reporting a problem

If you believe any material in this repository is incorrectly attributed or
should not be redistributed, please contact wcjohns@sandia.gov or
interspec@sandia.gov and it will be addressed.
