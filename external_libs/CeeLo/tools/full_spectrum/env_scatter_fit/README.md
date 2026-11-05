# Environmental-scatter fit campaign

This package implements the fit-first environmental-scatter workflow. It keeps
two deliberately different artifacts:

1. `CEEFIT-BLOCK-v1` worker records are transport handoff files. A chunked
   `CEEFIT-HDF5-v1` scientific archive keeps per-history sufficient
   statistics, 32 logical replicas, stopping covariance, component tags, exact
   lines, and the canonical configuration. It embeds the complete normalized
   raw records of each node's first four blocks; every later block is folded
   into the node's running sums as it is ingested and keeps a ledger row with
   the checksums its raw record would carry. Once an embedded record or fold
   and its materialized node have been durably verified, the original worker
   file is redundant and may be removed.
2. `CEEFIT-PAYLOAD-v1` is the dependency-free, checksummed runtime model. It
   contains only spline coefficients, domains, support thresholds, and locked
   validation metadata.

The worker format stores raw incident fluence. Fits are made to
`s_hat = raw_fluence / direct_fluence`; runtime evaluation requires the caller's
analytic `direct_fluence_cm2 = T(E0)/(4*pi*d_cm^2)` and returns both `s_hat` and
incident-fluence maps. Detector response is never part of either physics
artifact.

## Dependencies

Campaign design, allocation, block validation, registered axes, and rebinning
use only Python 3.10+ standard-library modules.

```bash
# Positive spline fitting
python -m pip install -r tools/full_spectrum/env_scatter_fit/requirements-fit.txt

# Fitting plus HDF5 archive ingestion
python -m pip install -r tools/full_spectrum/env_scatter_fit/requirements-archive.txt
```

Missing optional dependencies produce an actionable error; importing the core
package does not import NumPy or h5py.

The command can also be installed directly from the repository without an
internet-connected build environment:

```bash
python -m pip install --no-build-isolation tools/full_spectrum/env_scatter_fit
# Add [fit] or [archive] when dependency resolution is available.
ceefit --help
```

## Campaign workflow

Run every command from the CeeLo root (the directory that contains `tools/`)
so `tools.full_spectrum.env_scatter_fit` is importable. The C++ bake programs
live beside this package in `tools/full_spectrum/` and are built into
`build-rel/tools/full_spectrum/`, which is where the default worker command and
the examples below expect them:

```bash
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release
cmake --build build-rel --target env_kernel_bake env_fit_payload_eval \
  env_registered_e0_gate env_fit_block_tool env_kernel_merge env_scatter_export
```

Run an entire bake -- transport, ingestion, fitting, and evidence -- from one
build of one source tree. Code versions are not recorded; only the canonical
configuration identity (`config_sha256`) binds the artifacts together.

```bash
# Deterministic 10,000-point digitally shifted Sobol pool; forced G4,
# threshold, near/far, and domain-face anchors; disjoint 256/64/64 roles.
python -m tools.full_spectrum.env_scatter_fit design \
  --mode field --material soil \
  --output soil_campaign.json \
  --executable build-rel/tools/full_spectrum/env_kernel_bake

# One persistent process per node runs its initial four 262,144-history blocks.
python -m tools.full_spectrum.env_scatter_fit commands soil_campaign.json \
  --output run_initial_blocks.sh

# `--ungrouped` instead emits one independently schedulable process per block.

# Check a block before moving it from a worker.
python -m tools.full_spectrum.env_scatter_fit inspect-block ceefit_blocks/train/.../block_1.ceefit

# Validate and merge blocks. Duplicate IDs, history ranges, RNG seeds,
# configurations, or axes are rejected before mutation.
python -m tools.full_spectrum.env_scatter_fit ingest --output soil_fit.h5 \
  ceefit_blocks/*/*/*.ceefit

# Fit every supported continuum component and exact-line angular model by role,
# select rank/level on validation nodes, and write a candidate, validation
# report, and common-currency adaptive utility report.
python -m tools.full_spectrum.env_scatter_fit fit-archive \
  --archive soil_fit.h5 --campaign soil_campaign.json \
  --output soil.ceefitp --report soil_acceptance.json \
  --acquisition-utilities utilities.json

# Record the fit result. A failed validation report does not itself choose more
# histories or nodes; it advances the round into a unified-allocation state.
python -m tools.full_spectrum.env_scatter_fit campaign-round \
  --manifest soil_campaign.json \
  --validation-report soil_acceptance.json --archive soil_fit.h5 \
  --candidate soil.ceefitp.candidate \
  --output soil_campaign_round_1.json

# Read-only sizing: the blocks and archive growth that would bring every node
# to its noise target (the same assessment, summed over all nodes rather than
# one pass). Writes only the optional JSON report.
python -m tools.full_spectrum.env_scatter_fit project-continuation \
  --manifest soil_campaign_round_1.json --archive soil_fit.h5 \
  --adaptive-utilities utilities.json --output soil_projection.json

# Compare a complete 64-node/four-block acquisition batch with up to 64
# continuations (each at most doubling its node) in the same predicted
# validation-loss/CPU units.
# The utility report is byte-bound to the fit report, archive, model, and round.
python -m tools.full_spectrum.env_scatter_fit continue-archive \
  --manifest soil_campaign_round_1.json --archive soil_fit.h5 \
  --adaptive-utilities utilities.json \
  --output soil_campaign_allocated.json \
  --commands-output run_adaptive_work.sh

# A pass can only come from the strict fit-archive report, independently
# checked against the exact archive and serialized candidate. Two reports for
# consecutive round identities stop acquisition and lock the campaign state.
python -m tools.full_spectrum.env_scatter_fit campaign-round \
  --manifest soil_campaign_round_1.json \
  --validation-report soil_acceptance.json --archive soil_fit.h5 \
  --candidate soil.ceefitp.candidate \
  --output soil_campaign_round_2.json

# Generate the identity-bound detector comparison from the transport-complete
# pilot nodes, using conservative adjacent-pair aggregation for the 32-bin side.
python -m tools.full_spectrum.env_scatter_fit mu-fold-report \
  --archive soil_fit.h5 \
  --hpge-response responses/hpge.npz --nai-response responses/nai.npz \
  --output soil_mu_folds.json

# Decide whether the 64-bin pilot permits a 32-bin production mu axis. Without
# the generated identity-bound HPGe/NaI comparison, the decision is retain_64.
python -m tools.full_spectrum.env_scatter_fit mu-decision \
  --archive soil_fit.h5 --detector-validation soil_mu_folds.json \
  --output mu_decision.json

# The normal first pass has no external evidence yet and writes a candidate.
# Exercise that exact serialized candidate through both implementations. The
# deterministic bank covers physical faces, interior points, and values just
# below/at/above both the brems and pair-production thresholds.
python -m tools.full_spectrum.env_scatter_fit parity-report \
  --candidate soil.ceefitp.candidate \
  --cpp-evaluator build-rel/tools/full_spectrum/env_fit_payload_eval \
  --output evidence/python_cpp_parity.json

# Fold that same candidate and the archive's exact current 64-node final set.
# The response NPZ contract below keeps exact atoms off the continuum grid.
python -m tools.full_spectrum.env_scatter_fit detector-fold-report \
  --candidate soil.ceefitp.candidate --archive soil_fit.h5 \
  --campaign soil_campaign.json \
  --hpge-response responses/hpge.npz --nai-response responses/nai.npz \
  --output soil_detector_folds.json

# If detector-fold residuals cannot yet be judged because replica uncertainty
# exceeds 0.25%, continue only the named final nodes by one disjoint block per
# scheduling pass. Do not retire them while this is only an MC-noise failure.
python -m tools.full_spectrum.env_scatter_fit continue-archive \
  --manifest soil_campaign.json --archive soil_fit.h5 \
  --detector-validation soil_detector_folds.json --roles final \
  --output soil_final_continuation.json \
  --commands-output run_final_continuation.sh

# Calibrate all stored uncertainty estimators directly from immutable block
# records. Formal evidence needs at least sixteen blocks (eight independent
# adjacent block-pair comparisons); blocks from additional nodes add
# comparisons without reusing histories.
python -m tools.full_spectrum.env_scatter_fit calibrate-uncertainty \
  --archive soil_fit.h5 \
  --output evidence/uncertainty_calibration.json

# Compute the GEANT4 check from the archive's raw scalar N/S1/S2 and at least
# two independently seeded Livermore sphere-flux CSV arms per anchor. Each
# manifest arm binds relative CSV, GDML, and macro paths plus their SHA-256.
# Formal sphere-flux mode itself rejects non-Livermore EM, invalid seeds, and
# production cuts outside (0,0.7] mm. The output also records all
# geometry/scorer/run inputs.
python -m tools.full_spectrum.env_scatter_fit geant4-anchor-report \
  --archive soil_fit.h5 --manifest evidence/geant4_anchor_manifest.json \
  --output evidence/geant4_anchors.json

# Run the full-coverage registered-E0 gate through its integrated generator.
# The result records the exact argv, input parameters, raw 12x8 node maps, and
# output semantics. A failed report remains useful diagnostic evidence but
# cannot satisfy publication.
python -m tools.full_spectrum.env_scatter_fit registered-e0-report \
  --campaign soil_campaign.json \
  --executable build-rel/tools/full_spectrum/env_registered_e0_gate \
  --output evidence/registered_e0.json --threads 8

# Assemble the four detailed reports. Their paths are stored relative to the
# aggregate, and their exact file bytes are re-hashed again during promotion.
python -m tools.full_spectrum.env_scatter_fit assemble-evidence \
  --candidate soil.ceefitp.candidate --archive soil_fit.h5 \
  --campaign soil_campaign.json \
  --registered-e0-report evidence/registered_e0.json \
  --uncertainty-calibration-report evidence/uncertainty_calibration.json \
  --geant4-anchor-report evidence/geant4_anchors.json \
  --python-cpp-parity-report evidence/python_cpp_parity.json \
  --output evidence/soil_publication_evidence.json

# After folds/evidence arrive, promote that exact model without rerunning NMF.
python -m tools.full_spectrum.env_scatter_fit promote-candidate \
  --candidate soil.ceefitp.candidate --archive soil_fit.h5 \
  --campaign soil_campaign.json --detector-validation soil_detector_folds.json \
  --publication-evidence soil_publication_evidence.json \
  --output soil.ceefitp --report soil_promotion.json

# If an adequately resolved populated final set fails a physics gate, retire
# those IDs and replace them with 64 disjoint nodes from the unused pool.
python -m tools.full_spectrum.env_scatter_fit refresh-final \
  --manifest soil_campaign.json --failure-report soil_acceptance.json \
  --output soil_campaign_fresh_final.json
```

### Bounded-disk local workflow

`run-streaming` is the preferred path for a long local campaign. It runs a
persistent per-node worker group on the internal SSD while one background
writer compresses the preceding group into the canonical archive. The spool
limit is a soft, group-granularity bound: one worker group may temporarily
exceed it, but queued unarchived groups cannot grow without bound.

```bash
# The archive can live on the external drive; keep the short-lived spool on SSD.
python -m tools.full_spectrum.env_scatter_fit run-streaming \
  --manifest soil_campaign.json \
  --archive /Volumes/Backup/CeeLoWork/soil_fit.h5 \
  --spool scratch/ceefit_spool/soil \
  --spool-max-gib 4 --delete-packed-blocks
```

For each packed batch, the archive is closed and fsynced, every newly embedded
raw record is re-hashed, every newly folded block's ledger row is matched to
its worker block, and every affected node aggregate is checked exactly against
its fold state and raw records. Only then are the redundant `.ceefit` sources
unlinked, and the fold generation the batch replaced deleted. A failed worker, write, unplug, checksum, or closure check leaves the
source blocks in the spool. Re-running the same manifest uses a small,
stat-bound `<archive>.inventory.json` sidecar to skip archive-resident block IDs
without walking thousands of compressed HDF5 groups. A missing or stale index
is rebuilt from the archive. Partial spool output resumes without changing
histories or RNG seeds. Archive/spool location, thread scheduling, and shard
layout therefore do not change the Monte Carlo result.

Raw records and fold states use gzip level 1 without shuffle (materialized
node arrays keep shuffle). On the field/soil baseline's records that is 8%
smaller than with shuffle and faster to write; shuffle spreads the many exact
zeros of the sparse replica maps across byte planes. The setting changes
storage only, not Monte Carlo histories, precision, or fitting inputs. Archives
are created with persistent free-space tracking, so the space a replaced fold
generation or a folded raw record frees is reused instead of growing the file.

### Folded blocks

A node's materialized row is a fold over its blocks in ascending block-ID
order: compensated (Kahan) float64 sums of every moment, replica, stopping and
path-channel tensor, exact integer sums of the counts, and fsum/max/join merges
of the per-block scalars. `/folds/<node>/<k>` keeps that fold's complete state
after blocks 1..k -- every (total, compensation) pair, every integer sum, and a
ledger row per block (checksums, history and seed ranges, scalars) with a hash
chain over the block checksums -- sealed by `state_sha256`. Resuming the fold
from the state reproduces the fold over all blocks bit for bit, so:

- blocks 1-4 of every node (`RETAINED_RAW_BLOCKS`) stay raw records, for the
  uncertainty calibration;
- every later block that directly extends a node's fold is folded as the
  packer ingests it and is never written as a raw record; a block beyond a gap
  in the node's block sequence stays raw until the gap is filled;
- each ingest writes a new generation and moves the node's `current` pointer
  once the generation is complete; generations it replaced are deleted only
  after a post-write verification passes, and a stored state is checked
  against its seal whenever it is loaded, so an interrupted or failed ingest
  never costs the last verified fold and a damaged state is never extended;
- the archive snapshot, the materialized rows, the continuation schedule, the
  resume inventory and replay protection are the same as with every block
  raw (raw records and ledger rows feed one block inventory);
- the seal re-verifies a fold's checksum and ledger chain and recomputes the
  node's row from it; it cannot re-derive a folded block's contribution from
  that block's raw bytes, which were never stored.

A node's fold state measured at most 4.3 of its raw records at the 488-block
history cap (2 at four blocks), folding pseudo-blocks resampled from each
node's baseline; real blocks that fill cells the baseline never hit can push a
state toward its uncompressed size (~140 MB per node in production).
`project-continuation` prints both: the expected growth at five raw records per
node, and the ceiling at the uncompressed size, which the bake's round disk
guard enforces. `ingest --keep-raw-blocks` (`fold_blocks=False`) stores every
block raw, for tests and diagnostics.

By default the command performs one complete raw/materialized closure audit at
the end and writes `<archive>.seal.json`. Use `--no-seal` for several back-to-
back acquisition batches, then seal once before fitting or generating evidence:

```bash
python -m tools.full_spectrum.env_scatter_fit seal-archive \
  --archive /Volumes/Backup/CeeLoWork/soil_fit.h5
python -m tools.full_spectrum.env_scatter_fit verify-archive-seal \
  --archive /Volumes/Backup/CeeLoWork/soil_fit.h5
```

The seal records the archive scientific snapshot, canonical input JSON and
configuration hash, node/block inventory, and local file identity. The optional
`--full-file-hash` adds a bytewise HDF5 checksum; it is not needed for normal
local stale-result protection and costs another complete sequential read.

The closure audit costs about 0.3 s per raw record (single-threaded) plus a
fold-state checksum per folded node. Every seal
also writes `<archive>.verified-nodes.json`: per node, a closure identity
hashing its block inventory (block ids and record checksums, from raw records
and fold-ledger rows), its current fold generation and that generation's
seal, and its materialized row. `run-streaming --incremental-seal` (and `seal-archive
--incremental`) re-verifies only the nodes whose identity changed since that
ledger, runs every archive-wide check, and assigns the same snapshot identity;
the seal then says `"verification": "incremental-..."`. Adaptive rounds use it.
Publication (`promote-candidate`) accepts only a full seal and otherwise
re-verifies everything itself; the bake driver gives the archive a full seal
before publishing.

Fitting does not need `/blocks` or `/folds`. `fit-archive` reads the sealed
canonical archive directly (the bake driver does). When the archive lives on a
slow external disk, make a disposable aggregate-only cache on the internal SSD
and fit from it instead; both give the same model:

```bash
python -m tools.full_spectrum.env_scatter_fit fit-cache \
  --archive /Volumes/Backup/CeeLoWork/soil_fit.h5 \
  --output scratch/soil_fit_cache.h5
python -m tools.full_spectrum.env_scatter_fit fit-archive \
  --archive scratch/soil_fit_cache.h5 --campaign soil_campaign.json \
  --output soil.ceefitp --report soil_acceptance.json \
  --acquisition-utilities utilities.json
```

The cache copies the exact component/total moments, all 32 replica maps,
stopping moments, axes, diagnostics, and manifest; it changes no binning,
precision, weights, fit inputs, or acceptance calculation. It omits only the
embedded block records and is bound by a sidecar to the sealed canonical
snapshot. It can be deleted and recreated at any time. A measured production
soil archive needs about 2.9 GiB for this cache rather than about 13 GiB for the
canonical HDF5 or about 57 GiB for redundant worker files. Calibration, raw
closure, and publication evidence must still use the canonical archive, never
the fit cache.

The campaign manifest records every logical block and replica range/seed even
when four blocks are grouped into one persistent worker process. Continuation
uses the next per-node block sequence number. After a failed bound validation
round, `continue-archive --adaptive-utilities` is the production decision point.
It compares the summed predicted validation-loss reduction divided by measured
CPU for a complete 64-node/four-block acquisition batch with the same quantity
for up to 64 continuations. An existing node receives at most as many blocks as
it already has per scheduling pass, and never more than its noise projection
asks for, so a noisy pilot cannot commit a long speculative run while a node
that needs N times its histories reaches them in about log2(N) passes; each
continuation's reduction and CPU are priced for the blocks it receives.
Incomplete four-block baselines and validation-truth MC noise are mandatory
continuations: model refinement is not adjudicated against unresolved truth.
Joint-shape history stopping is evaluated after exact adjacent-cell summation
onto a predefined 12 registered-energy x 8 mu sentinel grid (or 1x8 for an
exact line). This is an allocation resolution, not a data reduction: native
384x64 masses and replicas remain in the archive, drive the fit, and retain
their native final-TV and detector-fold gates. Reports include native replica
TV as a diagnostic. On the pilot, cells below `1e-5` probability contribute
only 0.13 percentage points to a 6.77% median native noise floor, so
mass-conserving sentinel rebinning is used rather than thresholding nearly
empty cells.

Exact-line stopping is contribution-aware.  The environmental elastic-E0 atom
is combined with the deterministic direct E0 atom for allocation, final-node,
and detector-fold acceptance, while its standalone scattered contribution is
still archived and fitted.  The annihilation-511 atom always has an absolute
SE/error bound of 0.1% of the corresponding total.  Standalone 0.25% MC-noise
and 1% fit-error bounds activate only when its conservative three-sigma upper
contribution is at least 1% of that total.  Thus a weak atom cannot consume the
campaign through an irrelevant relative-error gate, but statistical
uncertainty cannot be used to classify a potentially material atom as weak.

The built-in new-node predictor uses normalized four-dimensional `q`, a
Matérn-3/2 localization of measured per-validation-node model loss, a novelty
(conditional-variance) penalty near existing training nodes, and an
eight-neighbor CPU-time model. It writes `CEEFIT-ADAPTIVE-UTILITY-v1`, including
new-node reductions and each existing node's localized validation loss. The
continuation branch combines that loss with the archive-derived MC-noise
reduction from one more block. This scheduling surrogate never changes
transport weights or physics.
The continuation evidence records the measured model-residual contribution and
the mandatory MC-precision guard contribution separately. They are summed only
for the common loss-reduction-per-CPU branch decision, so a sub-threshold model
residual is not relabelled as observed model error.

If any selection or final gate fails, `fit-archive` writes
`soil.ceefitp.candidate`, records the exact corrective policy, and exits nonzero;
it does not overwrite the requested production artifact. Failed final IDs are
retired into diagnostics/training only after archive and detector replica noise
are adequate enough to establish a physics failure. An uncertainty-only failure
instead schedules the next disjoint blocks at the same final IDs.
`refresh-final` draws 64 replacements from the persisted, previously unused
Sobol reserve only for an adjudicated physics failure. HPGe and NaI folds require
external detector-response matrices. Their report records every final node and
the response-matrix hashes and is bound to the configuration, candidate model,
and sorted raw-block archive snapshot. `CEEFIT-PUBLICATION-EVIDENCE-v1` further
requires the full-coverage registered-E0 gate, repeated-null S2 and
covariance calibration, GEANT4 anchor agreement, and Python/C++ parity. The
registered-E0 evidence must use schema `CEEFIT-REGISTERED-E0-GATE-v1`, carry
the same canonical `config_sha256`, and demonstrate all 12x8 registered joint
cells (including zero-estimate cells), total relative SE <=0.25%, and absolute
SE <=0.1% of total for significant joint cells and energy/angular marginals.
Its acceptance fields are `max_integrated_relative_error`, `max_joint_tv`,
`max_energy_marginal_tv`, and `max_angular_marginal_tv`; the
`registered-e0-report` command above generates it.
The campaign preregisters at most three total final generations: the initial
set and two fresh replacements. Exhausting that allowance requires a new
campaign revision instead of repeatedly testing until a favorable set passes.
Missing or stale evidence can only produce a candidate. Production loading rejects
candidate payloads unless the caller explicitly opts into diagnostic loading.

Publication also applies a structural production-configuration gate to the
archive and campaign themselves. The production design is exactly one
10,000-point candidate pool with disjoint initial roles of 256 training, 64
model-selection validation, and 64 untouched final nodes; every initial node
starts with four 262,144-history blocks. The archive requires the canonical
384-bin component-native axes, an exact 10-keV cutoff, the complete origin/order
and exact-atom inventories, and every baseline transport switch (next-event,
Rayleigh, bremsstrahlung, annihilation, attenuation) enabled. Its locked
transport values are a 150-keV brems electron floor, Doppler broadening off,
`max_vertices=50`, `roulette_tau=8`, `weight_floor=1e-5`, and a 4-cm scoring
radius. The radius must also close to the independently hashed GEANT4 geometry.
A 64-bin angular archive is accepted directly; 32 bins require the
identity-bound passing pilot decision. Reduced designs, fixtures, or altered
physics remain useful diagnostic fit inputs, but neither `fit-archive` nor
`promote-candidate` can publish them.

The untouched final 64 are the maximin-farthest points left after train and
validation selection. They are intentionally a sparse-region/worst-case stress
set, not an unbiased estimate of average-domain generalization. History caps are
reachable only in complete blocks: 244 blocks = 63,963,136 histories for
training/diagnostic nodes and 488 blocks = 127,926,272 for validation/final
nodes. A canonical 384x64 block is about 39.85 MB, so the initial 384-node,
four-block campaign is roughly 61 GB before HDF5 compression; plan local disk
and worker return bandwidth accordingly.

Fitting is deliberately bounded-memory. Continuum moment arrays are read one
physical component at a time, the 32 float32 replica maps are reduced in
eight-node row batches, and validation/final truth is scored one node at a
time. Archives chunk the cell dimension in at most 8,192-cell pieces. This
changes only I/O and memory scheduling: stored sufficient statistics, replica
sums, NMF inputs, spline candidates, and acceptance calculations remain the
same. Do not run several independent fit processes concurrently on a 16-GB
machine.

Payload size is diagnostic metadata, not a scientific acceptance gate. An
advisory 30 MiB value is reported only as a comparison. The local HDF5 Monte
Carlo archive has no payload-size constraint and is never distributed; the runtime product is the fitted functional form, so its natural
size is governed by the selected rank/spline complexity and validated accuracy.

Every formal block, archive, campaign, fitted candidate/payload, and generated
evidence report records the exact configuration identity (`config_sha256`).
Code versions are not recorded, so a whole bake must come from one build of one
source tree. Fitted candidates additionally record, under `fit_parameters`,
canonical input parameters, fit settings, output observable/normalization/edge
semantics, and per-model rank, level, and coefficient count.

The 64-to-32 pilot is similarly fail-closed: total joint-shape change and every
per-node HPGe/NaI total, spectral-TV, and exact-line fold change must be at most
0.2%, with response and archive hashes recorded. Otherwise production retains
64 angular bins.

External report contracts are intentionally explicit:

- `CEEFIT-MODEL-SELECTION-VALIDATION-v1` is emitted by `fit-archive` and binds
  every component/line and MC-noise gate to `config_sha256`,
  `archive_snapshot_sha256`, `model_content_sha256`, the exact eligible
  validation-node IDs, the next campaign round, and the current final-set
  generation. `campaign-round` reopens the archive and candidate before a pass
  can be recorded. Its canonical report hash and all binding identities are
  retained in `validation_history`; a report cannot be replayed as the next
  round. The only report-free route is the explicitly named
  `--diagnostic-validation-failed`, which always records a failure and can
  never increment the consecutive-pass counter.
  Every report and candidate also carries the exact
  `CEEFIT-CONTRIBUTION-AWARE-MODEL-SELECTION-v1` policy. Continuum-component
  relative gates activate at a conservative three-sigma contribution of 5%,
  exact-line relative gates at 1%, and weaker terms retain absolute
  0.1%-of-total yield and shape-L1 bounds. Elastic-E0 selection includes the
  analytic direct atom. `campaign-round` and `promote-candidate` reject reports
  or payloads with different gate semantics.
- `CEEFIT-DETECTOR-FOLD-VALIDATION-v1` carries `config_sha256`,
  `archive_snapshot_sha256`, `model_content_sha256`, the exact 64
  `final_node_ids`, and HPGe/NaI records. Each detector record has a
  `response_sha256`, canonical geometry/provenance hashes, explicit generator
  and response semantics, response dimensions, a closed conservative
  discretization-bound summary, and one metrics record per final node.
- `CEEFIT-MU-RESOLUTION-VALIDATION-v1` binds the 64-to-32 comparison to the
  configuration, archive snapshot, exact transport-complete pilot IDs, and
  response hashes plus the same geometry/provenance/generator evidence.
  `mu-fold-report` derives every continuum, spectral, exact-line, and
  discretization metric; a hand-authored assertion is not accepted.
  The embedded production decision retains absolute paths to the pilot archive
  and detector-fold report; campaign creation and production validation reopen
  and re-hash both sources, so the decision checksum is not merely format-level.
- `CEEFIT-PUBLICATION-EVIDENCE-v1` binds its four mandatory checks to that same
  configuration/archive, the model-content hash, and the fresh final IDs. Each
  check includes a relative path and SHA-256 of its detailed source report.
  Assembly computes that digest from the exact report bytes; validation opens
  the source again and rejects a missing, moved, edited, identity-mismatched, or
  metric-divergent report. A plausible-looking 64-character digest in a
  hand-written aggregate is not accepted.
- `CEEFIT-PYTHON-CPP-PARITY-v1` binds the serialized candidate model hash,
  configuration, and candidate-file hash, plus the exact request/result hashes,
  deterministic query inventory, compared-value count, and maximum absolute
  discrepancy. It compares continuum totals, every additive component, both
  incident and dimensionless masses, and all active exact-line angular masses.
- `CEEFIT-UNCERTAINTY-CALIBRATION-v1` binds the canonical configuration and
  sorted raw-block archive snapshot. Its top-level publication fields are
  `repeated_null_runs`, `s2_passed`, `coarse_covariance_passed`,
  `replica_covariance_passed`, and `family_coverage_passed`; all are derived by
  the command rather than supplied as operator assertions. `s2_passed`
  includes the scalar total and
  every stored continuum-component, total-continuum, line-component, and
  total-line cell, not only an integrated-yield check.
- `CEEFIT-GEANT4-ANCHOR-AGREEMENT-v1` is generated from raw, source-hashed
  48x32 sphere-flux CSV arms and the archive's own scalar moments. It requires
  `G4EmLivermorePhysics`, at least two byte-distinct independent arms per
  anchor with distinct RNG seeds, unbiased isotropic 4pi GPS sampling, exact
  GDML/macro hashes, matching energy/geometry/material/scoring-radius metadata,
  the archive energy cutoff (10 keV by default), and a 60 m or larger
  half-world for the infinite-field comparison. Every arm must contain raw
  `scattered_total_N/S1/S2` formed by summing the whole event's scattered sphere
  track length before squaring. The reference error is the larger of the
  merged per-history total-moment SE and the between-arm SE; cellwise errors
  are never added in quadrature for an integrated observable. Capped or
  transport-failed archive nodes cannot pass. A hand-authored
  `anchor_count`/`passed` pair cannot satisfy this contract.

  Formal inputs come from `env_scatter_export` and the
  `ceelo_g4val --sphere-flux` harness. Sphere-flux files without the exact
  total moments and GDML/macro identity metadata remain diagnostics and are
  rejected as publication evidence.

### Uncertainty calibration contract

`calibrate-uncertainty` reads only immutable `/blocks/<node>/<block>` records,
not rematerialized node aggregates or fold states. Numeric block IDs are paired
once within each node, so at least sixteen suitable blocks are needed for the
mandatory eight independent repeated-null comparisons. A folded node is
calibrated on its retained raw blocks 1-4 only (any other raw record it holds is
transient), so each node contributes two comparisons and a campaign's
calibration evidence is its baseline's (768 comparisons for the production
design). It tests total-fluence raw
per-history `S2` with standardized chi-square residuals and checks the marginal
`S2` of every cell in all four moment tensors. Physically supported continuum
components and exact-line segments must each be stochastic in at least eight
independent comparisons; unsupported/empty cells remain explicit deterministic
checks and cannot silently disappear from coverage. The exact 16-vector
stopping outer moments use covariance-rank-aware Mahalanobis residuals.
Cell coverage is accumulated in bounded memory, and its extreme-residual gate
uses an explicit family-size correction rather than a fixed cutoff that would
almost surely reject a correct production archive containing millions of
residuals.

Calibration blocks must collectively populate every physically required
origin/order and exact-line family. In particular, the pilot must provide eight
independent block-pair comparisons at or above the 1021.9979-keV pair threshold:
sixteen raw blocks at one node (an archive ingested with `--keep-raw-blocks`),
or the same number of disjoint blocks distributed among multiple
supra-threshold nodes -- in a folded archive, at least four such nodes. A node with no observed annihilation
continuum or 511-keV atom is retained as a valid failed coverage experiment,
not treated as proof of that family's uncertainty. Cell S2 arithmetic is
relative to its own operands, carries standard errors directly to avoid
subnormal squaring underflow, and uses only the IEEE square-underflow limit for
otherwise unrepresentable contributions. Float32 replica-map closure uses its
documented storage-precision bound, which remains far below the MC targets.

All four 32-replica map families are calibrated with finite-sample two-sample
Hotelling statistics. Deterministic projections include component and line
yields, twelve predefined total registered-energy bands, component coarse
rebinnings, exact-line Legendre moments, and smooth nonnegative response-like
folds. The report stores per-tensor and per-segment coverage plus one hash over
the complete projection inventory.

The strict JSON report records 68%, 95%, and 99.7% coverage diagnostics with
Wilson intervals, every paired block ID, unused odd blocks, omnibus
probabilities, all 16 stopping dimensions, and the deterministic projection
hash. Inadequate replication or family coverage, a malformed or singular
stochastic estimate, a residual outside covariance support, or any explicit
probability/residual gate writes a valid failed report and makes the command
exit nonzero. Pairing, the report, and its archive snapshot hash are invariant
to HDF5 group creation and block merge order.

### Detector response NPZ contract

`detector-fold-report` reads each response once as bytes, hashes those exact
bytes with SHA-256, and then opens the in-memory bytes with
`numpy.load(..., allow_pickle=False)`. A response is a standard uncompressed or
compressed NumPy NPZ containing exactly these arrays; extra keys are rejected:

| NPZ key | dtype and shape | Contract |
|---|---|---|
| `schema` | scalar string | `CEEFIT-DETECTOR-RESPONSE-v1` |
| `detector_name` | scalar string | exactly `HPGe` or `NaI`, matching the command-line slot |
| `response_semantics` | scalar string | `probability-per-incident-cell-mass-v1` |
| `line_interpolation` | scalar string | `linear-exact-energy-v1` |
| `geometry_json` | canonical scalar JSON string | exact `CEEFIT-DETECTOR-GEOMETRY-v1` detector/material/dimensions/orientation/source-surface record |
| `provenance_json` | canonical scalar JSON string | exact `CEEFIT-DETECTOR-RESPONSE-PROVENANCE-v1` transport, physics, histories, RNG, timestamp, and input-hash record |
| `generator_name` | scalar string | nonempty response-generator identity |
| `generator_version` | scalar string | nonempty response-generator version |
| `generator_sha256` | scalar string | lowercase SHA-256 of the exact generator artifact |
| `continuum_energy_edges_keV` | float64 `[Ne+1]` | strictly increasing incident-energy edges covering every final-node continuum axis |
| `mu_edges` | float64 `[Nmu+1]` | strictly increasing incident-direction edges spanning exactly `[-1,1]` |
| `spectrum_edges_keV` | float64 `[Ns+1]` | strictly increasing recorded-spectrum edges |
| `continuum_response` | float64 `[Ne,Nmu,Ns]` | nonnegative probability in each recorded bin per incident continuum cell mass |
| `line_energy_grid_keV` | float64 `[Nl]` | strictly increasing exact incident-line energies covering every supported final-node atom |
| `line_response` | float64 `[Nl,Nmu,Ns]` | nonnegative angular response at each exact-line grid energy |

The last axis of each response tensor must sum to at most one for every
incident cell. HPGe geometry must name active material `HPGe`; NaI geometry
must name `NaI(Tl)`. Both use
`cosine-to-detector-axis;mu=1=head-on-v1`, explicitly binding response angle to
the archive coordinate. `dimensions_cm` is a nonempty positive-valued object.
Provenance requires nonempty transport/version/physics/RNG information,
positive histories per response grid point, and at least one bound input
SHA-256.

Continuum archive masses and all 32 replica maps are
mass-conservingly rebinned to the response incident axes before folding. The
candidate is evaluated directly on those same edges. For an exact E0 or
511-keV atom, `line_response` is linearly interpolated at the atom's exact
energy and folded against its angular masses; no continuum bin centre is used.

When response boundaries split an archive cell, both report generators compute
a distribution-free response-kernel-diameter bound: each cell's unknown
sub-cell allocation is allowed to be any convex combination of every response
kernel it overlaps. Yield, normalized-spectrum, and exact-line metrics are
inflated by that bound. The bound itself must be at most 0.2% for the 64-to-32
decision and 1% for untouched-final validation. Thus a fine-looking fold cannot
pass by relying on an undocumented uniform-within-cell approximation. Exact
incident-axis alignment makes the corresponding bound identically zero.

For every one of the exact 64 fresh final IDs and for both detectors, the
report records continuum folded-total relative error, normalized spectral TV,
exact-line angular-fold relative L1 error, and a conservative replica scalar:
the largest detector-folded yield relative SE or spectral-shape replica-TV SE
across continuum and every supported atom. Structural coverage, dimensions,
support, history accounting, finite and nonnegative data, candidate/archive
identity, or line-energy failures abort
without replacing an existing report. The report separates `physics_passed`
from `noise_passed`: excessive replica uncertainty causes same-node
continuation, while a physics failure retires the final set only after the
0.25% detector-fold noise target has been met.

`mu-fold-report` applies the same response contract to every eligible 64-bin
pilot node. Its fine side uses the archived 64-bin continuum and exact-line
angular masses. Its candidate side first sums adjacent pairs exactly into 32
bin-integrated masses, then conservatively maps those masses to the response
axis. The discretization interval covers response boundaries splitting either
the 64-bin reference cells or the wider 32-bin candidate cells. It atomically
emits folded-total relative change, normalized spectral TV, and exact-line
angular-fold relative L1 change for HPGe and NaI. Those are the three metrics
consumed by `mu-decision`; missing coverage, provenance, or a bound above 0.2%
keeps the production campaign at 64 bins.

## Data and model contracts

- Node IDs are SHA-256 over `CEEFIT-NODE-v1\0`, the ASCII configuration hash,
  and `(E0_keV,d_cm,hsrc_cm,hdet_cm,extent_cm)` as five little-endian binary64
  values. The leading eight digest bytes are interpreted big-endian.
- Block IDs are positive per-node sequences. A block has 262,144 histories in
  32 equal logical replicas; replica seeds hash the complete
  `(config,node,block,replica)` tuple. The formal Python archive reader rejects
  every other block size, even when the binary checksum and replica arithmetic
  are otherwise valid.
- HDF5 has required `/manifest`, `/axes`, `/nodes`, `/continuum`, `/lines`, and
  `/diagnostics` groups plus immutable per-block records (`/blocks`) and
  per-node fold states (`/folds`). Scientific arrays use gzip-compressed
  chunks. Aggregate node arrays are rebuilt from the node's fold state and the
  raw records beyond it, in numeric block-ID order after every ingestion, so
  continuation is bitwise independent of block arrival order and of how many
  blocks were folded. Elastic-E0 and annihilation-511 form a fixed two-line
  axis; below pair threshold the absent 511 line has `N=histories`, zero
  `S1/S2`, zero replicas, and a false support bit.
- The reader recomputes canonical configuration, node identity, replica ranges
  and seeds, registered/physical axes, component taxonomy and support, line
  inventory, diagnostic channels, moment inequalities, replica closure, and
  stopping-vector closure before archive mutation. A valid outer checksum is
  not sufficient to bypass these semantic checks.
- The first 12 stopping entries use the 0.1%-of-total absolute-SE target only
  when their band carries at least 1% of total fluence. The four signed
  Legendre moments and replica joint-shape sentinels use 0.25%-of-total; tiny
  cells never receive a relative-error requirement. Every raw block and HDF5
  `/stopping` group carries the same validated `CEEFIT-STOPPING-v1` canonical
  JSON definition: exact registered-axis index ranges for bands 0–11 and the
  standard Legendre orders 1–4, including angular and normalization semantics.
- Exact-line allocation and acceptance use the same contribution-aware policy:
  combined analytic-direct plus environmental elastic E0, and an always-active
  0.1%-of-total absolute bound for annihilation-511 with standalone relative
  yield/shape gates activated by its three-sigma 1% materiality bound.
- Component fits use positive log-yields and softmax geometry weights on a
  four-dimensional clamped tensor cubic basis. Joint energy-angle archetypes
  are nonnegative mixtures of individually unit-integral cubic B-spline PDFs.
  NMF is performed on cell-integrated masses, then projected through the exact
  archive-cell integration design using nonnegative multiplicative NNLS.
- Runtime integration evaluates those spline PDFs over requested cell edges;
  archive bins are never treated as delta functions or merely renamed splines.
- Every payload declares its complete continuum/line `model_inventory`, the
  `brems_activation_threshold_keV`, an `accepted_for_publication` bit, and a
  `model_content_sha256`. The content hash covers the executable model and its
  physics/domain metadata but intentionally excludes later validation and
  publication state, allowing an external detector report to bind the exact
  candidate before that same payload is promoted.
- Every model carries the full deterministic candidate-pool `q` box and rejects
  extrapolation. The payload separately enforces rectangular physical bounds
  and the coupled `|hsrc-hdet| < distance` constraint.

## Performance policy

Production workers keep every validated transport process enabled. The accepted
exact refactors are a persistent per-node worker pool, reusable geometry-trace
scratch, collision cross-section reuse, and one detector-leg material trace
shared by the Compton and Rayleigh attenuation evaluations. The isolated
200,000-trace benchmark measured the shared trace as bit-identical with a 3.16x
speedup; campaign timing remains part of the recorded node diagnostics.

Primary antipodal/hemisphere pairing and near-field detector-cone antithetic
sampling are implemented only behind benchmark modes. Their bounded tests do
not meet the paired-bias and variance-times-wall-time gates, so both stay off.
None of this alters Rayleigh, bremsstrahlung, next-event estimation,
attenuation, exact atoms, or the normal source-direction sampler. Analytic
single scatter is not a baseline replacement and remains only a future,
closure-gated control variate candidate.

## Tests

Run from the CeeLo root. Tests that need NumPy or h5py are skipped when the
corresponding extra is not installed, and the tests that execute
`env_fit_payload_eval` or `env_registered_e0_gate` are skipped until those
binaries exist in `build-rel/tools/full_spectrum/`.

```bash
python -m unittest discover -s tools/full_spectrum/env_scatter_fit/tests -t .
```

Tests include cross-language node/seed vectors, CEEFIT block corruption checks,
checksum-valid semantic-corruption rejection, incremental/permuted HDF5
merge-order invariance, folded archives reproducing every-block-raw archives
bit for bit (rows and snapshot), gap handling, fold-state corruption and
missing-record detection, verification before deleting a replaced fold,
compressed chunk checks, mass-conserving rebinning,
uncertainty gates and repeated-null covariance calibration, deterministic
adaptive design/fresh-final selection, genuine spline integration, and
byte-for-byte runtime payload fixtures shared with the C++ evaluator.
