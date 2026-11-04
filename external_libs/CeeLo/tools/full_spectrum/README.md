# Full-spectrum bake

Everything needed to regenerate, verify and accept CeeLo's four full-spectrum
payload families from this source tree alone:

| Family | Driver | GEANT4 | Final files (in `data/full_spectrum/`) |
|---|---|---|---|
| HPGe detector response | `bake/detector.sh` | no | `<detector>_detector_response/` (one per detector) |
| Source escape | `bake/source_escape.sh` | optional parts | `source_escape/` |
| Charged-particle source | `bake/charged_particle.sh` | no | `charged_particle_source/` |
| Environment scatter | `bake/env_scatter.sh` | no | `env_scatter_<tier>/` |

A *bake* runs one family's driver in production mode: it generates the Monte
Carlo truth banks, fits the payload, scores it against held-out banks and
independent checks (the *gates*), and copies the accepted payload into
`data/full_spectrum/<family>/`. Nothing a bake produces is committed except
the final payload of a bake whose gates all passed (section 6).

Run every command from the CeeLo root (the directory holding `CMakeLists.txt`
and `tools/`).
Every driver resumes where it stopped: after an interruption or a fixed error,
re-run the same command.

Contents:
1. [Set up the machine](#1-set-up-the-machine)
2. [Build and verify](#2-build-and-verify)
3. [Smoke run](#3-smoke-run)
4. [Production bakes](#4-production-bakes)
5. [Re-running, resuming and changing a recipe](#5-re-running-resuming-and-changing-a-recipe)
6. [Committing the results](#6-committing-the-results)
7. [Troubleshooting](#7-troubleshooting)
8. [What is in this directory](#8-what-is-in-this-directory)

## 1. Set up the machine

A bake machine needs:

- a C++17 compiler, CMake 3.20 or newer, Eigen 3 and Boost (unit test framework);
- Python 3.10 or newer with the packages in `requirements.txt` (numpy, scipy,
  h5py, matplotlib);
- GEANT4 11.4.0 with its datasets (about 2 GB), only for the optional
  GEANT4 parts of the source-escape bake (section 4.2);
- disk and memory for the families you bake:

| Family | Disk (work area) | Memory |
|---|---|---|
| Detector | about 0.05 GB | 2 GB |
| Source escape | about 1 GB | 4 GB |
| Charged particle | about 1.5 GB | 8 GB |
| Environment scatter, per tier | 15 GB archive after the baseline, about 25 GB after the adaptive rounds (at most ~70 GB), + 4 GB spool | 16 GB |

Bake each family on one machine with one build: different machines, compilers
or build flags give different Monte Carlo streams, and an environment-scatter
archive must come from a single build.

Keep the work area (`WORK`, and the environment's `ARCHIVE_DIR`) on an APFS
volume, not ExFAT or FAT. Programs built or installed on a Mac carry the
`com.apple.provenance` extended attribute and pass it to every file they
write; on a non-native volume macOS stores it in a `._<name>` AppleDouble file
beside each one, which the bake's directory scans would read as data. Copy
finished results to such a disk with `cp -X`. `bake/preflight.sh` checks this.

The instructions are for macOS (Apple silicon or Intel) with Homebrew; the
same steps work on Linux with the distribution's packages.

```bash
brew install cmake ninja eigen boost python@3.12
python3.12 -m venv ~/ceelo-venv
~/ceelo-venv/bin/pip install -r tools/full_spectrum/requirements.txt
export PY=~/ceelo-venv/bin/python      # every driver uses $PY
```

(On Ubuntu 24.04: `sudo apt-get install build-essential cmake ninja-build
libeigen3-dev libboost-test-dev python3-numpy python3-scipy python3-h5py
python3-matplotlib`, and `PY` can stay `python3`.)

GEANT4 11.4.0, if it is not installed yet (built from source; about an hour):

```bash
curl -fsSL https://github.com/Geant4/geant4/archive/refs/tags/v11.4.0.tar.gz | tar -xz
cmake -S geant4-11.4.0 -B geant4-build -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$HOME/geant4-11.4.0 -DGEANT4_BUILD_MULTITHREADED=ON \
    -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_GDML=OFF -DGEANT4_USE_QT=OFF
cmake --build geant4-build && cmake --install geant4-build
```

Then, in every shell that builds or bakes:

```bash
source $HOME/geant4-11.4.0/bin/geant4.sh
export PY=~/ceelo-venv/bin/python
cd <InterSpec checkout>/external_libs/CeeLo
```

Long bakes run for hours: start them with `nohup ... &` (or inside `tmux` or
`screen`), keep the Mac from sleeping with `caffeinate -is`, and give them the
lowest CPU priority with `nice -n 20` so the machine stays usable:

```bash
nohup caffeinate -is nice -n 20 tools/full_spectrum/bake/<driver>.sh > <driver>.out 2>&1 &
```

Section 4.0 runs every family in one such command.

## 2. Build and verify

```bash
tools/full_spectrum/bake/build.sh
```

It checks the prerequisites, configures and builds a Release tree in
`build-rel/` (or `$BUILD`), runs the C++ unit tests (`ctest`), the Python tests
of the environment-scatter package and of the bake tools
(`tools/full_spectrum/tests`), and the C++/Python lockstep checks. It
takes about 8 minutes on an M1 Pro (compiling 1 minute, C++ tests 3,
Python tests 3). It ends with `BUILD_DONE`. **Stop if anything fails.**
The GEANT4 programs are built only when `geant4.sh` was sourced before the
first configure; they are needed only for the source-escape bake's optional
GEANT4 parts. To add them later: source `geant4.sh`, delete
`build-rel/CMakeCache.txt` and run `build.sh` again.

Dependencies CMake cannot find by itself go in `CMAKE_ARGS`, for example
`CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/opt/deps" tools/full_spectrum/bake/build.sh`.

Then measure the machine (about 5 minutes; run it on an otherwise idle machine):

```bash
tools/full_spectrum/bake/timing_probe.sh
```

Divide each time by the reference below to get the machine's speed factor,
and scale the production times in section 4 by it:

| Probe | Reference time |
|---|---|
| detector cell 661.657 keV, 500k events, 1 thread | 3.0 s |
| detector cell 2614.511 keV, 3M events, 1 thread | 115 s |
| source-escape cell 3000 keV Z=26 tau=4, 1M events, 1 thread | 11 s |
| electron cell 5000 keV Z=82 1 g/cm2, 60k events, 1 thread | 0.8 s |
| env node 2614 keV at 5 m, 4 x 262,144 histories, all cores | 15 s |
| env node 59.5 keV at 5 m, 4 x 262,144 histories, all cores | 5.0 s |

(Reference: Apple M1 Pro, 10 cores, with other work running at the same time,
so an idle M1 Pro is somewhat faster. The two detector cells run with spectrum
fidelity; their times are the measured fidelity-off references scaled by the
on/off CPU-time ratio, 0.80 and 2.29. The production times in section 4 are
estimates built from these per-cell costs.)

## 3. Smoke run

A smoke run executes every step of every driver on tiny banks. It checks the
machine, the build and the scripts end to end before any long run. Its payloads
are meaningless and its gates are expected to fail (they are reported, not
enforced). It writes only under `bake_work/smoke/`.

```bash
SMOKE=1 tools/full_spectrum/bake/detector.sh
SMOKE=1 tools/full_spectrum/bake/source_escape.sh
SMOKE=1 tools/full_spectrum/bake/charged_particle.sh
SMOKE=1 tools/full_spectrum/bake/env_scatter.sh
```

Each must end with its `..._DONE` line (`DETECTOR_BAKE_DONE`,
`SOURCE_ESCAPE_BAKE_DONE`, `CHARGED_PARTICLE_BAKE_DONE`,
`ENV_SCATTER_ROUND_DONE`) and exit status 0. Smoke times on an M1 Pro: detector 30 s,
source escape 2.5 min, charged particle 20 s, environment scatter 13 min.
Afterwards `rm -rf bake_work/smoke`.

## 4. Production bakes

The four families are independent; bake them in any order, one at a time (each
uses every core). Each driver prints a timestamped line per step and keeps
every step's full output in `bake_work/<family>/logs/<step>.log`. When a gate
fails in production the driver stops with the gate's log; the gate table in it
says which bound failed.

### 4.0 Everything, unattended

On a machine left to bake (the environment campaign's adaptive rounds alone
run for days), one command does all of it:

```bash
nohup caffeinate -is nice -n 20 env ENV_ROUNDS=40 ARCHIVE_DIR=/Volumes/Big/ceelo_env \
    tools/full_spectrum/bake/all.sh > bake_all.out 2>&1 &
tools/full_spectrum/bake/status.sh      # any time, from any shell
```

`all.sh` runs `preflight.sh` (work-area file system and free space, memory,
Python, power, priority; a problem stops it), `build.sh` (a failure stops it),
the timing probe once, then each detector in `DETECTORS` (default
`"hpge_coax detective_x"`), source escape, charged particle, and environment
scatter with `ROUNDS=$ENV_ROUNDS`. A family that fails does not stop the
others; each driver's output is in `$WORK/<stage>.out`, each stage's exit
status and last `_DONE`/`ERROR` line in `$WORK/all/stages.txt`, and the closing
`status.sh` report in `$WORK/bake_report.txt`. `FAMILIES` selects families.
Re-run the same command to resume: each driver skips the steps it has done.
`SMOKE=1 tools/full_spectrum/bake/all.sh` is the whole smoke run of section 3
(with `SMOKE=1` the work area defaults to `bake_work/smoke`).

`ENV_ROUNDS=0` (the default) stops the environment campaign at its baseline
checkpoint (section 4.4); the campaign ends itself after two consecutive
passing rounds, so a full campaign gives a number it will not reach. Each
round checks its disk space before it transports anything (section 4.4).

Settings shared by every driver (environment variables):

| Variable | Default | Meaning |
|---|---|---|
| `JOBS` | cores - 1 | concurrent Monte Carlo cells |
| `BUILD` | `build-rel` | build tree (section 2) |
| `PY` | `python3` | Python interpreter |
| `WORK` | `bake_work` | banks, fits, reports, logs |
| `RESULTS` | `data/full_spectrum` | where accepted payloads are copied |
| `REDO` | (empty) | steps to run again although done (section 5) |
| `SMOKE` | 0 | 1 = smoke run (section 3) |
| `PRECISION` | `final` | `fast` / `medium`: 1/8 / 1/2 of the events, for a driver with tiers (section 4.1) |

### 4.1 HPGe detector response

```bash
tools/full_spectrum/bake/detector.sh                        # hpge_coax
DETECTOR=detective_x tools/full_spectrum/bake/detector.sh   # any other detector
```

A detector payload describes one detector; bake each one you need. `DETECTOR`
names it (default `hpge_coax`):

- `hpge_coax`, the campaign's built-in reference: a germanium crystal 8.0 cm
  in diameter and 8.0 cm long with a 0.5 cm-radius, 5.5 cm-deep bore, dead
  layers of 0.05 / 0.08 / 0.08 cm (front / side / back) and a 1 mm aluminium
  endcap.
- any other name `N` reads `tools/full_spectrum/bake/config/N.geometry`, a
  `CEELO_GEOMETRY 1` file (section 8, "Detector geometry"). To add a detector,
  write that file -- from an InterSpec DRF that carries a `<CeeLoGeometry>`
  block, `$PY tools/full_spectrum/drf_geometry.py DRF.xml --output
  tools/full_spectrum/bake/config/N.geometry` -- and run the driver with
  `DETECTOR=N`. `detective_x.geometry` is a placeholder until the real
  dimensions are known: the repository's ORTEC Detective-X-like crystal
  (6.5 cm diameter, 5.0 cm long, 0.45 cm x 3.8 cm bore, 0.07 cm dead layers)
  inside the Micro-Detective's housing layers.

The bank manifests record the detector's label and geometry text, the fitter
embeds both in the payload, and the evaluator derives the rate-surface anchor,
the feature map's passive-layer landmark (55 keV for `hpge_coax`, 62.2 keV for
`detective_x`) and its geometric acceptance from them. Each detector has its
own work area (`bake_work/detector_N/`), seeds and published directory. The
driver:

1. Places the bank energies (`detector_bank_energies.py`). Energies below
   100 keV are written for `hpge_coax` and scaled by E_L / 55 keV, so the
   low-energy banks sit at the same place in the passive-layer turn-on of any
   detector (E_L from `ceelo_rate_surface_anchor --info`; `hpge_coax`'s
   55.000006 keV keeps the written energies exactly, every energy being
   rounded to 0.05 keV). A held-out energy (audit and arrival banks) within
   2% of any training energy is moved by the smallest multiple of 0.5% that
   clears them all; for `hpge_coax` that moves audit_a's 1000 / 2000 / 5000 /
   10000 keV to 980 / 2040 / 5100 / 10200 and audit_b's 70 / 120 keV to
   68.6 / 119.4. The placed recipe is written to `bake_work/detector_N/recipe.txt`;
   a work area baked with another recipe is refused (with the diff) unless
   `REDO` names the steps to run again (section 5).
2. Generates nine truth banks with `ceelo_detector_response_campaign`
   (seeds derive from the detector, the bank's salt and the cell), with
   detector bremsstrahlung, fluorescence and spectrum fidelity on
   (`--detector-spectrum-fidelity`: Compton K-vacancy x-rays and the Detailed
   electron walk, as `EfficiencyCalculator::enable_spectrum_fidelity`). A
   history whose deposit is the energy less an escaped crystal K x-ray is a
   K-escape atom whatever made the vacancy, so those lines reach the payload
   at every energy, not only below ~300 keV. A
   fixed-angle cell is a parallel beam onto the detector's envelope at one
   energy and incidence angle; an arrival cell (`--incidence CLASS`) is an
   isotropic field over a band of direction cosines -- iso [-1, 1], front
   [0, 1], back [-1, 0], side [-0.3, 0.3] -- each photon drawing its own
   direction and its origin on the bounding disk, both redrawn on a miss, so
   it is the area-weighted average of the fixed-angle response: what a
   scattered field is folded with. Cells run until their statistics target
   between the minimum and maximum event counts.

   | Bank | Energies (keV) | cos(theta) or class | Cells | Events per cell | Role |
   |---|---|---|---|---|---|
   | training | 19 from 10 to 12000 | -1 to 1, step 0.25 | 171 | 0.2-0.5 M | fit |
   | dense_energy | 10 between the training energies | 8 | 80 | 0.3-0.6 M | fit |
   | dense_angles | 11, 50 to 12000 | 16 | 176 | 0.25-0.6 M | fit |
   | highstat | 8 from 1332 to 12000 | 5 | 40 | 3 M | fit |
   | lowE | 8 from 0.64 to 1.32 E_L (passive-layer turn-on; 35-72.5 for `hpge_coax`) | 9 | 72 | 0.6-3 M | fit |
   | validation | 11 | 8 | 88 | 0.25-0.6 M | gates |
   | audit_a, audit_b | 10 each, off every training grid | 8 each, to grazing | 80 + 80 | 0.6-3 M | reported |
   | arrival | 15: 0.75-1.64 E_L and 110-850 | iso, front, back, side | 60 | 4 M | gates |

3. Fits the payload (`fit_detector_response.py`: one bandwidth triple, the
   pair-family and 511-keV yield fields widened or narrowed per component,
   and every yield field's nodes solved so its kernel average reproduces the
   training cells, `--yield-node-iterations 16`). The fit's payload, scored on
   the validation bank, is the candidate. The fitter refuses arrival cells
   (`cell_kind detector_arrival`), so they cannot leak into training. Then
   `ceelo_detector_finite_distance` measures the candidate's finite-distance
   depths (a point source on the axis at 5 cm to 10 m, 21 energies, 1M events
   per point, up to 16x where a point is noisy; 25-60 minutes on a 10-core Mac,
   most of it the 5-12 MeV points) and writes them into
   it: `predict`'s acceptance for a source at a finite distance.
4. Scores the candidate on both audit banks, and on the arrival bank
   (`score_arrival_bank.py`: the payload averaged over each class's band by
   `ceelo_detector_model_evaluate --incidence`, Gauss-Legendre in the polar
   angle per hemisphere, against the sampled truth).
5. Runs 4M-event direct Monte Carlo spot checks at normal incidence
   (122, 661.657 and 1332.492 keV) and evaluates the candidate there.
6. Gates (`check_detector_gates.py`). On the validation bank: detected yield
   error median <= 0.15% and p95 <= 0.5%; pulse-shape total variation median
   <= 1.3% and p95 <= 4%; full-energy atom error median <= 0.5% and p95 <= 1.5%;
   each spot-check full-energy yield within 0.5% of direct Monte Carlo. On the
   arrival bank, in every class and energy band (40-100, 100-400, 400-1000
   keV): detected yield error median <= 0.15% and max <= 0.5% beyond 2 sigma
   of the cell's own MC noise; continuum total variation (without the
   zero-deposit and full-energy atoms) median <= 3x its MC noise floor; valley
   density (the continuum between the Compton edge and the full-energy line)
   error median <= 3% and max <= 10%; the envelope's mean projected area
   within 5 sigma of the MC boundary's. The arrival full-energy atom is
   reported until its misfit is fixed.
7. Publishes `detector_response_payload.txt`, the candidate's `report.json`
   and `gates.json` to `data/full_spectrum/N_detector_response/`.

Time: the banks are 0.6-1.3 G events at 20-45 microseconds each, 4-14
CPU-hours (about 1-1.5 hours on a 10-core Mac; spectrum fidelity is ~2.7x of
it, nearly all above 2 MeV), and the fit and scoring add about
10 minutes. The banks are almost all of it; every later step takes
minutes, so a failed gate never means regenerating banks (section 5).

`PRECISION=fast` (1/8 of every bank's and spot check's events) or `medium`
(1/2) is a real bake at lower precision, to see where a detector stands in a
fraction of the time. It writes only under `bake_work/precision_<tier>/`,
reports its gates without enforcing them (their bounds assume final noise
floors) and publishes into that work area's `results/`. Only the detector
driver has the tiers; the others refuse to run below `final`, and `all.sh`
then defaults `FAMILIES` to `detector`.

**Known state:** no detector passes every gate yet. `hpge_coax` misses two:
the validation full-energy p95 (2.25% against 1.5%, from 11 MeV) and the arrival
front 40-100 keV detected-yield median (0.163% against 0.15%); `detective_x_fit`
also misses the arrival valley medians (3.3-3.4% against 3%). See TODO.md, "The
detector-response full-energy atom between nodes" and "The detector bakes miss
arrival gates". The banks are still worth generating: fix the fitter, then
`REDO="fit candidate finite_distance audit_a audit_b arrival spot" tools/full_spectrum/bake/detector.sh`.

### 4.2 Source escape

```bash
tools/full_spectrum/bake/source_escape.sh
```

The payload describes the photons escaping a point source at the centre of a
shielding sphere (energy, effective Z, areal density, hydrogen fraction). It
is assembled from four fitted parts, each from its own truth banks (all
`ceelo_source_escape_campaign`, Rayleigh and Compton Doppler on, shield K and L
x-rays and electron-impact ionization on):

| Part | Banks | Fit |
|---|---|---|
| Secondaries (pair annihilation, electron bremsstrahlung, shield K-alpha / K-beta / L x-rays) | 1122 cells (11 energies x 17 Z x 6 tau, 1M histories); 75 off-grid tau cells (4M); 672 hydrogen-mixture cells (1M) | `fit_source_secondaries.py`, after `check_kernel_span.py` accepts the mechanistic bases |
| Elastic line | the secondary and hydrogen banks | `fit_elastic_escape.py` (one scale, Compton damping pinned to 1) |
| Three-or-more Compton scatters | 5880 training/validation cells (42 energies x 20 Z x 7 tau, 100k histories, low-count cells refined to 2M) and 792 test cells (200k), each with its mechanistic-model spectrum | `fit_multiple_scatter_surrogate.py` |
| Photonuclear (optional) | GEANT4: 152 cells at 8-12 MeV plus 11 off-grid cells (4M each) | `fit_photonuclear.py` |

`assemble` writes `source_escape_payload.txt` with `ceelo_source_model_evaluate
--write-payload` and checks that the written file reproduces the assembled
model exactly. Then the payload is scored on data no fit selected on:

- the off-grid tau bank (`score_offgrid_holdout.py`), including the
  fluorescence atoms and the share of Compton-scattered shield x-rays, which
  the model does not produce;
- the deep-tau spurious tail (`check_deep_tau_tail.py`);
- the 150-cell acceptance grid, {59, 100, 250, 661, 1000, 3000} keV x
  Z {6, 13, 26, 48, 82} x {0, 1, 5, 15, 50} g/cm2, with a 20M-history CeeLo
  truth per cell (`run_source_anchor.py`, `score_acceptance_grid.py`); a cell
  is scored when its collided leakage is at least 5e-5 per primary, which
  leaves out the unshielded and the opaque cells whatever the statistics.

`check_kernel_span.py` gates the mechanistic bases before any fit. Its bars
apply to the hydrogen-free rays; the hydrogen-mixture rays are gated
separately, with their own bounds for the two components whose bases are
weakest there.

The gates (`check_source_escape_gates.py`) bound each score by the last
accepted payload's value with headroom, max(1.25 x last, last + 0.5 points),
unless a documented bar exists; `gates.json` records each bound's origin. They
also re-run the C++/Python lockstep checks against the build.

Two parts are optional GEANT4 runs (`source <geant4>/bin/geant4.sh` first):

- `PHOTONUCLEAR=1` bakes the photonuclear term. It matters only above about
  8 MeV and evaluation leaves it off by default, so a default bake omits it
  from the payload. It costs about 450 CPU-hours (1-4 ms per GEANT4 event at
  8-12 MeV, 650M events): about two days on a 10-core Mac.
- `GEANT4_CHECK=1` adds the acceptance grid's GEANT4 reference leg (100k
  histories per cell, a few CPU-hours) and gates the payload against it too.

Either can be added to a finished bake later by re-running with the setting:
the driver re-runs exactly what it affects.

Time, estimated: about 12 CPU-hours for a default bake (multiple-scatter
banks about 7, including a 1-4 s model evaluation per cell; secondary banks
about 2; acceptance grid about 1; fits 1-2): 2-3 hours on a 10-core Mac. Disk: about 1 GB, most of it the multiple-scatter banks.

### 4.3 Charged-particle source

```bash
tools/full_spectrum/bake/charged_particle.sh
```

The payload describes the photons escaping a monoenergetic electron or
positron at the centre of a sphere (kinetic energy T, effective Z, thickness
in units of the particle's CSDA range r). Truth is `ceelo_source_electron_campaign`
(512 bins, 10 keV photon floor; shield K and L x-rays and electron-impact
ionization on, the unscattered x-rays tallied by family as yield-only atoms whose
bases come from `ceelo_charged_particle_model_evaluate --fluorescence-bases`):

| Bank | Grid | Events per cell | Cells |
|---|---|---|---|
| training, sub-grid A (both particles) | T 50-12000 keV (13) x r 0.02-10 (11) x 18 Z | 60k (Z 32, 70: 120k) | 5148 |
| training, sub-grid B (electrons) | T 120-10000 keV (18) x r 0.2-6 (16) x 15 Z | 200k | 4320 |
| training, sub-grid C (both particles) | T 240-900 keV (8) x r 0.1-2 (7) x 15 Z | 300k | 1680 |
| off-grid validation (electrons) | T (12) x Z {32, 70} x r (8), all off the training grid | 300k | 192 |
| off-grid final test (both particles) | T (11) x Z {16, 47, 79} x r (7) | electrons 400k, positrons 1.6M | 462 |

Z 35 and 64 are whole-Z holdouts: the fits never see them. Their positron
cells in sub-grids A and C carry 16x their sub-grid's events (960k and 4.8M):
at 1x the in-flight annihilation yields sit at their own MC noise floor, and
the holdout planes are where the positron gates are scored. The driver checks
each bank's cell count.

1. `fit_structured`: `fit_charged_particle_source_structured.py` fits the
   yield fields and a starting shape model with the frozen per-component
   structure in `bake/config/charged_particle_structure.json`.
2. `joint_extract`, `joint_build`: `fit_charged_particle_joint_shape.py` fits
   the final shape blocks jointly (16 rounds, settings frozen in
   `bake/config/charged_particle_joint_settings.json`), keeping the yield
   fields byte for byte.
3. `joint_score`: scores the payload on the holdout Z, the off-grid validation
   bank and, once, the off-grid final test.
4. `parity`: the C++ evaluator against the Python predictor on probe cells.
5. Gates (`check_charged_particle_gates.py`): payload format invariants, C++
   vs Python parity within 0.5%, and per-component yield and shape errors on
   the holdout planes and the off-grid final test, bounded by the last
   accepted payload's values with headroom. A component is scored in a cell
   when its truth yield is at least 25 counts at 60k events (4.2e-4 photons per
   primary), whatever the cell's own statistics.

The two config files are hyperparameters chosen by selection studies; the bake
does not re-select them (`fit_charged_particle_joint_shape.py reselect`
re-derives the joint settings if the model changes).

Time: the banks are 2.99 G events at 13-30 microseconds each, about 1 hour
10 minutes on an M1 Pro (the 1.86 G-event recipe took 45 minutes there); the
fits add about 7 minutes. Disk: about 1 GB (about 30 KB per compressed cell,
plus the structured fit's 0.5 GB cache of parsed cells).

### 4.4 Environment scatter

Put the archive on a large disk, and keep the work area (spool and fit cache)
on a fast local disk:

```bash
ARCHIVE_DIR=/Volumes/BigDisk/ceelo_env tools/full_spectrum/bake/env_scatter.sh
```

`MODE=field MATERIAL=soil` is the default tier; the others are
`MATERIAL=concrete` and `MODE=air`. Each tier is its own campaign with its own
work area (`bake_work/env_<tier>/`).

The first run is the baseline:

1. `design`: a deterministic campaign of 256 training, 64 validation and 64
   final nodes (E0 20-3000 keV, distance 0.25-25 m, source height 0-4 m,
   detector height 0.3-4 m; 384 energy x 64 angle bins). The manifest
   (`campaign.json`, about 16 MB) binds every block ever transported for the
   campaign; the driver never regenerates it once it exists.
2. `streaming`: transports 1536 blocks of 262,144 histories each, packing them
   into `<ARCHIVE_DIR>/<tier>_archive.h5` (about 15 GB) through a 4 GB local
   spool, then verifies the archive's seal. It resumes after an interruption.
   These first four blocks of every node stay in the archive as raw records;
   every later block is folded into its node's running sums as it is packed
   (the archive keeps the fold's exact state and a ledger row per block), so
   the adaptive rounds add little to the archive.
3. `calibrate`: repeated-null calibration of the uncertainty estimators
   (`calibration.json` must say `"passed": true`).
4. `fit1`: every rank/level candidate, fitted from the sealed archive
   (`FIT_CACHE=1` fits from an aggregate-only copy in the work area instead,
   for an archive on a slow disk; the model is the same). Before the adaptive
   rounds converge, validation fails: `fit-archive` then exits nonzero and
   writes `<tier>.ceefitp.candidate`. That is expected.
5. `record1`: records round 1 in `round1/manifest.json`.

Time: about 2 hours for the baseline's transport on a 10-core Mac (the
single packing thread, about 19 s per node, is the bottleneck), plus 10-30
minutes for `fit-archive`. The driver then stops and prints what to review.

**Checkpoint before adaptive rounds:**

| Check | Where | Go if |
|---|---|---|
| Calibration | `calibration.json` | `"passed": true` for all four families |
| Model selection | `round1/acceptance.json`, `selection.<component>.max_yield_relative_error` and `.max_joint_tv` | trending toward yield error <= 0.5% and joint TV <= 1% |
| Noise projection | `project-continuation` (below) | most nodes projected under their history caps, and the archive growth fits the disk |
| Disk | `du -sh <ARCHIVE_DIR>` | about 15 GB, and more free than `project-continuation` prints |

The noise projection (read-only, about 30 s):

```bash
W=bake_work/env_field_soil
$PY -m tools.full_spectrum.env_scatter_fit project-continuation \
    --manifest $W/round1/manifest.json --archive <ARCHIVE_DIR>/field_soil_archive.h5 \
    --adaptive-utilities $W/round1/utilities.json
```

prints, per role, how many nodes already meet their noise target, how many need
their whole history cap, and the blocks and archive growth to bring every node
to its target. The field/soil baseline of the Oct 2026 bake projected +67,527
blocks (17.7 G histories), about 80 hours of transport at the baseline's rate
(63 hours of it in the adaptive rounds; the final nodes' share comes after the
campaign stops), and about +12.7 GB of archive, at most +55 GB: the folded
nodes' states, expected at about five raw records per node and at most their
uncompressed size. With every block kept raw it would have been about 445 GB.
An archive created before folding (it predates persistent free-space tracking)
never reuses the space of replaced fold states; `project-continuation` says so,
and a new baseline avoids it.

Then run adaptive rounds:

```bash
ROUNDS=3 ARCHIVE_DIR=/Volumes/BigDisk/ceelo_env tools/full_spectrum/bake/env_scatter.sh
```

Each round allocates more histories or nodes from the last round's utilities
(`continue-archive`; a continued node at most doubles its blocks per round),
transports them, refits and records the round. Before transporting, a disk
guard projects an upper bound on the round's archive growth from its scheduled
blocks (raw records for new nodes' first four blocks, fold-state growth for
continued nodes) and stops the driver unless `DISK_MARGIN_GIB` (default 20)
plus the spool stay free on the archive's and the work area's disks; free
space and re-run to resume. A round's seal re-verifies only the nodes it
changed. The
campaign stops itself after two consecutive passing rounds; the driver then
gives the archive a full seal (publication requires one), ends with
`ENV_SCATTER_BAKE_DONE` and publishes the candidate and its reports.

For the field/soil baseline that is about 27 rounds (simulated from the
round-1 projection: validation nodes first, 64 nodes per round), the largest
~6,400 blocks (~7.5 hours of transport). Packing, verification and sealing
handle each new block about once (rebuilding a node starts from its fold state,
not from all of its blocks): about 10 hours of single-threaded work over the
campaign.

**Known state:** turning a converged candidate into a published payload
(`promote-candidate`) also needs detector-fold, GEANT4-anchor and full-statistics
registered-E0 evidence, and the producers of the first two do not exist yet
(TODO.md).

## 5. Re-running, resuming and changing a recipe

- **Resume:** run the same command. Finished steps print `skip <step>`, and
  the Monte Carlo runners skip every finished cell.
- **Run a step again:** `REDO="<step> <step>"` (names as printed in the
  `run`/`skip` lines), or `REDO=all`. The source-escape driver re-runs every
  step downstream of a changed one by itself; the other drivers do not, so name
  those steps too. For example, after a fitter change:
  `REDO="fit candidate finite_distance audit_a audit_b arrival spot" tools/full_spectrum/bake/detector.sh`.
- **Start a family over:** delete `bake_work/<family>/` (for environment
  scatter also the archive).
- **What invalidates what:** banks depend on the transport physics (`src/`) and
  on the recipe in the driver; fits depend on the banks and the fitter. After
  a physics change, regenerate the banks (delete the family's work area). After
  a fitter change, re-run from the fit step.
- **Changing a recipe** (energies, event counts, fit settings) means editing
  the driver: the drivers are the record of how each payload was made. The
  detector driver records its placed bank recipe in the work area and stops
  when it differs, printing the diff; re-run with `REDO` naming the changed
  banks and every later step. A bank runner also refuses a bank directory
  holding cells its campaign does not plan (energies a recipe change moved):
  move those cells out, since the fitter and the scorers read every cell in
  the directory.

### Comparing against a previous payload

To see what a new bake changed, score the payload it replaces (`P` below; keep
a copy before the new bake publishes over it) on the new banks, with the same
tools the gates use, once the family's banks and fits exist. Then compare the
reports with the candidate's.

```bash
EVAL=build-rel/tools/full_spectrum/ceelo_source_model_evaluate
# Detector: the validation and audit banks (and score_arrival_bank.py
# --payload P on banks/arrival).
for bank in validation audit_a audit_b; do
  $PY tools/full_spectrum/fit_detector_response.py --evaluate-payload P \
    --validation bake_work/detector_hpge_coax/banks/$bank \
    --photoelectric-basis-size 16 --single-compton-basis-size 48 \
    --multiple-compton-basis-size 96 --pair-no-compton-basis-size 96 \
    --pair-single-compton-basis-size 64 --pair-multiple-compton-basis-size 96 \
    --output bake_work/detector_hpge_coax/previous_$bank
done
# Source escape: the off-grid holdout and the acceptance grid.
$PY tools/full_spectrum/score_offgrid_holdout.py --payload P --evaluator $EVAL \
  --offgrid bake_work/source_escape/banks/offgrid \
  --output bake_work/source_escape/previous_offgrid.json
$PY tools/full_spectrum/score_acceptance_grid.py --payload P --evaluator $EVAL \
  --grid bake_work/source_escape/banks/acceptance \
  --output bake_work/source_escape/previous_acceptance.json
# Charged particle: the holdout planes and both off-grid banks.
$PY tools/full_spectrum/fit_charged_particle_joint_shape.py score --payload P \
  --fit-dir bake_work/charged_particle/structured \
  --banks test,offgrid_validation,offgrid_final_test \
  --output bake_work/charged_particle/previous_report.json
```

## 6. Committing the results

Commit a family's payload only after its production driver ended with its
`..._DONE` line and every gate passed. Never commit smoke results, banks,
intermediate fits, or a candidate from an unfinished campaign.

```bash
git add -f data/full_spectrum/hpge_coax_detector_response
git commit
```

(`data/full_spectrum/` is git-ignored so nothing lands there by accident;
`-f` adds the published files deliberately.) The published `gates.json` is the
payload's acceptance record; mention the build machine and the bake dates in
the commit message.

## 7. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `missing .../tools/full_spectrum/<program>` | not built, or `BUILD` points elsewhere | run `bake/build.sh` (with the same `BUILD`) |
| `the GEANT4 environment is not set` | `geant4.sh` not sourced in this shell | `source <geant4>/bin/geant4.sh`, re-run |
| `GEANT4 datasets not found under ...` | datasets missing or not mounted | install them, or mount them (section 1b) |
| `missing .../g4_source_escape_campaign` with `PHOTONUCLEAR=1` or `GEANT4_CHECK=1` | built without GEANT4 | source `geant4.sh`, delete `<build>/CMakeCache.txt`, re-run `build.sh` |
| `missing Python modules: ...` | packages not installed for `$PY` | `$PY -m pip install -r tools/full_spectrum/requirements.txt` |
| `step <name> failed (exit N); its log is ...` | that step's program failed | read the end of the log; fix; re-run the same command |
| `gate <name> FAILED` | the payload misses an acceptance bound | the log's table shows which; see the family's section |
| exit status 137, or `Killed` | out of memory | close other applications, or bake on a machine with more memory |
| `unable to lock file` (HDF5) | file locking unsupported on the mount | `export HDF5_USE_FILE_LOCKING=FALSE` (one writer only ever opens the archive) |
| `run-streaming` very slow | spool on a slow or network disk | keep `WORK` on a fast local disk; only `ARCHIVE_DIR` goes on the big disk |
| env `fit-archive exit status 2` | validation not yet passing | expected before the adaptive rounds converge (section 4.4) |
| `._<name>` files in the work area, or `PREFLIGHT FAILED` naming exfat/msdos | work area on a non-native volume (section 1) | move `WORK`/`ARCHIVE_DIR` to APFS; copy results out with `cp -X` |
| `disk guard: the round would leave less than N GiB free` | the next environment round does not fit on the disk | free space or move `ARCHIVE_DIR` to a bigger disk, then re-run the same command |
| ExFAT or USB disk writes at tens of MB/s | the disk is shared with a Time Machine backup, or the FSKit ExFAT driver | keep `WORK` local; let the backup finish before the baseline's streaming |
| env `reusing the existing campaign manifest` | resuming a campaign | expected; the manifest is never regenerated |
| `runtime model inventory does not exactly match fitted models` | a component had no validation support (reduced design) | do not shrink the production design |


## 8. What is in this directory

**Drivers and settings** (`bake/`): `common.sh` (settings, the resumable
`step`/`gate`/`publish` helpers), `build.sh`, `timing_probe.sh`,
`detector.sh`, `source_escape.sh`, `charged_particle.sh`, `env_scatter.sh`;
`all.sh` (every family, unattended), `preflight.sh` and `status.sh`;
`bake/config/` holds the charged-particle fit's frozen hyperparameters and the
detectors' `<name>.geometry` files.
`requirements.txt` lists the Python packages.

**Monte Carlo truth** (C++, built into `<build>/tools/full_spectrum/`):

| Program | Produces |
|---|---|
| `ceelo_detector_response_campaign` | one detector cell: parallel beam onto the detector's envelope, or (`--incidence`) an isotropic field over a band of arrival cosines |
| `ceelo_source_escape_campaign` | one source-escape cell: point source in a shielding sphere |
| `ceelo_source_electron_campaign` | one charged-particle cell: electron or positron in a sphere |
| `env_kernel_bake` | environment-scatter fit blocks (next-event fluence estimator) |
| `g4_source_escape_campaign`, `g4_photonuclear_xs` | GEANT4 source-escape cells and photonuclear cross sections (built when GEANT4 is found) |

**Campaign runners** (Python; each writes one bank, one gzip-compressed CSV per
cell, `cells/<id>_ceelo.csv.gz`, skipping finished cells; `bank_io.py` is the
one place that names, opens and identifies cells, and it also reads
uncompressed `.csv` cells):
`run_detector_response_campaign.py`, `run_source_secondary_campaign.py`
(secondary, off-grid and hydrogen banks), `validate_multiple_scatter.py` and
`refine_multiple_scatter_cells.py` (multiple-scatter banks),
`run_photonuclear_campaign.py`, `run_source_anchor.py` (acceptance grid),
`run_charged_particle_source_campaign.py`, and the environment-scatter package
`env_scatter_fit/` (`python -m tools.full_spectrum.env_scatter_fit`; see its
README for every subcommand).

**Fitters:** `fit_detector_response.py`, `fit_source_secondaries.py`,
`fit_elastic_escape.py`, `fit_multiple_scatter_surrogate.py`,
`fit_photonuclear.py`, `fit_charged_particle_source_structured.py`,
`fit_charged_particle_joint_shape.py`; the environment-scatter fit is
`env_scatter_fit fit-archive`. Shared libraries: `fit_charged_particle_source.py`,
`compare_reference.py`, `validate_single_scatter.py`, `source_escape_groups.py`.

**Evaluators** (C++): `ceelo_detector_model_evaluate` (`--incidence`
averages the payload over an arrival band),
`ceelo_source_model_evaluate` (also assembles the source-escape payload),
`ceelo_charged_particle_model_evaluate`, `env_fit_payload_eval`, and the
helpers `ceelo_source_cross_section_features` (cross sections and features the
fitters share with C++) and `ceelo_rate_surface_anchor`.

**Checks and gates:** `check_detector_gates.py`, `check_source_escape_gates.py`,
`check_charged_particle_gates.py`, `check_kernel_span.py`,
`check_deep_tau_tail.py`, `score_arrival_bank.py`, `score_offgrid_holdout.py`,
`score_acceptance_grid.py`,
and the C++/Python lockstep guards `check_basis_lockstep.py`,
`check_feature_lockstep.py`, `check_base_lockstep.py`.

**Bank placement:** `detector_bank_energies.py` scales a detector bake's
sub-100 keV energies to the detector's passive-layer landmark and keeps the
held-out banks 2% off the training grid.

**Detector geometry:** `ceelo_detector_response_campaign --geometry FILE`
runs any detector described by a `CEELO_GEOMETRY 1` text file
(`GeometryDescriptor::to_text`); `drf_geometry.py` writes one from an
InterSpec DRF's `<CeeLoGeometry>` block, and `--print-geometry` writes a
built-in detector's. `ceelo_rate_surface_anchor --geometry FILE --info` prints
what the detector payload derives from it (landmark, K-line element and edge,
envelope).

**Not used by the drivers:** `env_kernel_merge` (merges environment-scatter
shards baked on several machines), `env_fit_block_tool` (inspects fit blocks),
`env_scatter_export` (GEANT4 inputs for environment-scatter validation),
`env_registered_e0_gate` (the registered-energy publication gate), and the
evaluator benchmarks `bench_detector_response_model` and
`bench_charged_particle_source`.
