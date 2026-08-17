# Pu 610–775 keV known-truth fixtures

This directory contains a deliberately small set of plutonium spectra for
relative-activity regression tests.  `manifest.json` is the machine-readable
source of truth for provenance, isotope values, file hashes, and test
tolerances.

## Included data

- Three EC JRC `UPu-Gamma` ORTEC microDetective spectra: the replicate-zero,
  four-hour, 10 cm, Fe 16 mm/Cd 0 mm measurements of CBNM Pu93, Pu70, and
  Pu61.
- The two distinct long background measurements referenced by those three
  spectra.  Each N42 file was made by summing the 192 five-minute SPC files in
  its JRC background session with SpecUtils; the Pu61 and Pu93 spectra share
  the 2017-10-11 background, so it is stored only once.
- IAEA International Database of Reference Gamma Spectra (IDB-v2024-01)
  plutonium record 1461.  It is useful for parsing and robustness checks but
  has no paired background and is not an unconditional fit-accuracy gate.

The foreground SPE files and IDB file are byte-for-byte copies with their
original basenames.  The fixture selection itself is the only change to those
files.  See the manifest for their original relative paths and SHA-256 hashes.

## Truth convention

The certificate values are recorded at their stated reference date and their
uncertainties are absolute two-sigma weight-percent values.  The five Pu
isotopes are normalized independently to total remaining Pu; Am-241 is kept as
a separate Am-241/total-Pu mass ratio.

Acquisition-time values in the manifest were recomputed from the certificate
date to the start of each spectrum.  The calculation uses the half-lives,
atomic masses, and Pu-241-to-Am-241 branching ratio in the repository's
`data/sandia.decay.xml`: exponential depletion for each Pu isotope, the
two-member Bateman equation for Am-241 ingrowth/decay, and normalization only
after decay.  The constants are copied into the manifest so the oracle is
auditable and does not silently move when the nuclear database changes.
For IDB record 1461, the separately published acquisition-date composition is
also preserved and checked against this derivation to within 0.02 absolute
weight-percent point per component.

The manifest also records a model-derived **effective age** for each JRC case.
It is obtained by inverting the two-member Pu-241-to-Am-241 ingrowth equation
at acquisition, with the idealized initial condition of zero Am-241 at one
chemical separation.  This is only the conditional age used by the known-age
fit-accuracy diagnostic; it is not a certified separation date or independent
truth about the sample's history.  The equation, constants, inputs, units, and
full-precision results are all recorded in `manifest.json`, and the fit test
loads those values from the manifest.

The SPE `$DATE_MEA` fields contain timezone-free, whole-second wall-clock
timestamps.  For the JRC files they are seven hours later than the explicitly
UTC acquisition timestamps in the accompanying JRC metadata.  The manifest
therefore keeps UTC acquisition truth separate from each embedded SPE header
value and records absolute live/real-time comparison tolerances that account
for the respective source fields' rounding.

Do **not** use the JRC metadata element named
`decay_corrected_mass_fractions_wt_pct` as acquisition truth.  The locally
generated metadata currently repeats the certificate-date table under that
name rather than performing a decay correction.  Likewise, the isotopics
embedded in the IDB SPE remarks are informative source metadata; the manifest
contains the consistently recomputed oracle.

## Attribution and license

The JRC spectra and background inputs are from the European Commission Joint
Research Centre `UPu-Gamma` open-data collection.  The exact upstream
copyright notice is preserved in `COPYRIGHT-JRC.txt`.  The data are licensed
under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/);
reuse must credit the European Union/EC JRC and indicate the background
summation described above.

IDB record 1461 is from the IAEA Nuclear Data Section's *International
Database of Reference Gamma Spectra*, version IDB-v2024-01.  The record names
IAEA as data provider and EC-JRC Karlsruhe as laboratory.  Cite the
[dataset DOI](https://doi.org/10.61092/iaea.gaef-dhff) and the
[data descriptor](https://doi.org/10.1038/s41597-025-05455-1).  IDB-v2024-01
is distributed under CC BY 4.0.
