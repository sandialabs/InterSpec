#!/usr/bin/env python3
# CeeLo: a Monte Carlo photon-transport library for computing gamma-ray
# and X-ray detector efficiency - developed as part of InterSpec.
#
# Copyright 2026 National Technology & Engineering Solutions of Sandia, LLC
# (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.
# For questions contact William Johnson via email at wcjohns@sandia.gov, or
# alternative emails of interspec@sandia.gov.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

"""Joint shape fit of the charged-particle source payload.

fit_charged_particle_source_structured.py fits each continuum component's
shape model in four independent stages: (a) per-cell NNLS of the truth
spectrum against the 32-term basis, (b) count-weighted NMF of those shape
vectors into R archetypes, (c) per-cell NNLS of the mixture weights, and
(d) penalized least squares of the log-odds fields.  Each stage minimizes a
different surrogate, so their errors compound.  This tool keeps the structured
payload's yield fields and domains byte for byte and replaces (b)-(d) with one
alternating minimization of the DEPLOYED shape error: the weighted squared
error of the model's grouped spectrum against the grouped MC truth, which is
what the reported `shape_total_variation` measures (up to L1 versus L2):

    E(U, L) = sum_i c_i || A_i (softmax(B_i L)_i U) - t_i ||^2 + L^T Pen L

  A_i  (G x 32)  grouped component basis of cell i (columns unit-sum)
  t_i  (G,)      unit-sum grouped MC truth
  B_i  (M,)      tensor B-spline design row on the weight grid
  U    (R x 32)  archetypes, U >= 0, rows sum to 1
  L    (R-1 x M) log-odds field coefficients, gauge l_{R-1} = 0

Stage (a) never enters: the target is the MC spectrum itself, so the basis's
own representation error is inside the objective rather than frozen into the
target.  The fit starts from the four-stage solution, and the best iterate by
deployed error is kept, so it cannot end worse than its start.

U-step: E is quadratic in U for fixed weights; its normal equations are
assembled with one matrix product, unit-sum rows are imposed by a large
quadratic penalty and nonnegativity by NNLS on a Cholesky factor.  L-step:
L-BFGS-B with the analytic softmax gradient.  Rescaling an archetype row is an
exact gauge of the parameterization (the evaluator renormalizes the mixture,
and a constant shift of a log-odds field is exact because the B-spline basis
is a partition of unity), so the payload's unit-sum-rows invariant is free.

Subcommands, in bake order (FIT is the structured fitter's --output):
  extract   per-component shape datasets from FIT's parsed-bank caches
  build     joint-fit the shape blocks and write the payload
  score     score a payload on the whole-Z holdout planes and off-grid banks
  parity    compare the C++ evaluator with the Python reference evaluator
  reselect  optional, not part of the bake: re-run the staged selection of
            the joint-fit settings (fits the bank's train split, scores its
            validation split plus the off-grid validation bank)
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import scipy.linalg
import scipy.optimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_charged_particle_source as fit  # noqa: E402
import fit_charged_particle_source_structured as structured  # noqa: E402

BASIS_SIZE = fit.BASIS_SIZE
PAYLOAD_NAME = "charged_particle_source_payload.txt"
SHAPE_COMPONENTS = tuple(
    name for name in fit.MODEL_COMPONENTS
    if name not in structured.YIELD_ONLY_COMPONENTS)
# Groups of the reported shape metric (fit.grouped_metrics).
GROUPS = 32
# Bank label -> parsed-cell cache written by the structured fitter.
CELL_CACHES = {
    "bank": "bank_cells",
    "offgrid_validation": "offgrid_validation_cells",
    "offgrid_final_test": "offgrid_final_test_cells",
}
SETTING_KEYS = (
    "rank", "weight_grid", "weight_lambda", "loss", "count_weighting",
    "init_lambda")
LOSSES = ("l2", "cell_l1", "group_l1")


def log(message: str) -> None:
    print(message, flush=True)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# The structured fit's output directory
# ---------------------------------------------------------------------------

def fit_holdout_z(fit_dir: Path) -> set[int]:
    report = fit_dir / "report.json"
    if not report.exists():
        raise SystemExit(
            f"{report} is missing: run "
            "fit_charged_particle_source_structured.py --output "
            f"{fit_dir} first")
    value = json.loads(report.read_text()).get("holdout_z")
    if value is None:
        raise SystemExit(
            f"{report} records no holdout_z: re-run the structured fit")
    return {int(item) for item in value}


def fit_cells(fit_dir: Path, label: str, holdout_z: set[int]) -> list:
    """Cells of one bank with their scoring split: the training bank is
    split by the structured fit's rule, an off-grid bank is one split named
    after the bank."""
    cells = structured.read_cell_cache(
        fit_dir / "cache" / f"{CELL_CACHES[label]}.pkl")
    for cell in cells:
        cell.split = (
            fit.split_for(cell.identity, cell.z_value, holdout_z)
            if label == "bank" else label)
    return cells


# ---------------------------------------------------------------------------
# Shape datasets
# ---------------------------------------------------------------------------

class ShapeData:
    """Per-(cell, component) rows written by the extract subcommand."""

    FIELDS = (
        "identity", "split", "feature", "log_yield", "sigma_rel",
        "counts", "shape", "matrix", "truth")

    def __init__(self, arrays: dict):
        self.__dict__.update(arrays)

    @classmethod
    def load(cls, path: Path) -> "ShapeData":
        if not path.exists():
            raise SystemExit(
                f"{path} is missing: run the extract subcommand first")
        with np.load(path, allow_pickle=False) as blob:
            return cls({key: blob[key] for key in cls.FIELDS})

    @property
    def size(self) -> int:
        return self.feature.shape[0]

    def subset(self, mask: np.ndarray) -> "ShapeData":
        return ShapeData({key: getattr(self, key)[mask] for key in self.FIELDS})


def grouped_operator(size: int) -> np.ndarray:
    """(groups x size) summation operator matching fit.grouped_metrics."""
    usable = min(GROUPS, size)
    boundaries = np.linspace(0, size, usable + 1, dtype=int)
    operator = np.zeros((usable, size))
    for index in range(usable):
        operator[index, boundaries[index]:boundaries[index + 1]] = 1.0
    return operator


def extract_component(cells, name, papf, shape_cache, basis_cache) -> dict:
    """The deployed shape metric depends only on the model's 32-vector
    shape: prediction = total * (basis @ shape), and grouping then
    normalizing cancels `total`.  So each (cell, component) row needs only
    A = G @ basis (G x 32) and the unit-sum grouped truth t.  Selection of
    rows mirrors the structured fitter's (>= 25 expected counts)."""
    rows = {key: [] for key in ShapeData.FIELDS}
    operators: dict[int, np.ndarray] = {}
    for cell in cells:
        if not fit.component_applies(cell, name):
            continue
        truth, variance = fit.truth_for(cell, name)
        total = float(truth.sum())
        base = fit.component_base(cell, name)
        if total * cell.events < 25.0 or base <= 0.0:
            continue
        basis_key = (
            name, cell.energy_keV, cell.z_value, cell.areal_density,
            cell.edges.size, float(cell.edges[-1]), papf)
        basis = basis_cache.get(basis_key)
        if basis is None:
            basis = fit.component_basis(cell, name, papf)
            basis_cache[basis_key] = basis
        key = (cell.identity, name, papf)
        coefficients = shape_cache.get(key)
        if coefficients is None:
            coefficients = structured.robust_nnls(basis, truth / total)
            shape_cache[key] = coefficients
        coefficient_sum = float(np.sum(coefficients))
        if coefficient_sum <= 0.0:
            continue
        operator = operators.get(truth.size)
        if operator is None:
            operator = operators[truth.size] = grouped_operator(truth.size)
        matrix = operator @ basis
        grouped = operator @ truth
        if matrix.shape[0] < GROUPS:
            # Zero rows on both sides add nothing to the objective.
            pad = GROUPS - matrix.shape[0]
            matrix = np.vstack([matrix, np.zeros((pad, BASIS_SIZE))])
            grouped = np.concatenate([grouped, np.zeros(pad)])
        rows["identity"].append(cell.identity)
        rows["split"].append(cell.split)
        rows["feature"].append(fit.cell_feature(cell, name))
        rows["log_yield"].append(math.log(max(total / base, 1e-300)))
        rows["sigma_rel"].append(
            math.sqrt(max(float(np.sum(variance)), 0.0)) / total)
        rows["counts"].append(total * cell.events)
        rows["shape"].append(np.asarray(coefficients) / coefficient_sum)
        rows["matrix"].append(matrix)
        rows["truth"].append(grouped / max(grouped.sum(), 1e-300))
    return {
        "identity": np.asarray(rows["identity"], dtype=str),
        "split": np.asarray(rows["split"], dtype=str),
        "feature": np.asarray(rows["feature"], dtype=float).reshape(-1, 3),
        "log_yield": np.asarray(rows["log_yield"], dtype=float),
        "sigma_rel": np.asarray(rows["sigma_rel"], dtype=float),
        "counts": np.asarray(rows["counts"], dtype=float),
        "shape": np.asarray(
            rows["shape"], dtype=float).reshape(-1, BASIS_SIZE),
        "matrix": np.asarray(
            rows["matrix"], dtype=float).reshape(-1, GROUPS, BASIS_SIZE),
        "truth": np.asarray(rows["truth"], dtype=float).reshape(-1, GROUPS),
    }


# ---------------------------------------------------------------------------
# Metrics on the deployed shape
# ---------------------------------------------------------------------------

def deployed_tv(matrix: np.ndarray, truth: np.ndarray,
                shapes: np.ndarray) -> np.ndarray:
    """0.5 * L1 distance between the grouped model and grouped truth
    spectra, both unit-sum: the reported `shape_total_variation`."""
    predicted = np.einsum("igk,ik->ig", matrix, shapes)
    total = np.maximum(predicted.sum(axis=1, keepdims=True), 1e-300)
    return 0.5 * np.abs(predicted / total - truth).sum(axis=1)


def summarize(values) -> dict:
    finite = np.asarray([value for value in values if np.isfinite(value)])
    if finite.size == 0:
        return {"count": 0}
    return {
        "count": int(finite.size),
        "median": float(np.median(finite)),
        "p95": float(np.percentile(finite, 95)),
        "maximum": float(np.max(finite)),
    }


def selection_objective(summary: dict) -> float:
    if summary.get("count", 0) == 0:
        return math.inf
    return summary["median"] + 0.3 * summary["p95"]


# ---------------------------------------------------------------------------
# Shape model and alternating fitter
# ---------------------------------------------------------------------------

class ShapeModel:
    def __init__(self, domain, weight_grid, archetypes, log_odds):
        self.domain = np.asarray(domain)
        self.weight_grid = tuple(weight_grid)
        self.archetypes = np.asarray(archetypes)
        self.log_odds = None if log_odds is None else np.asarray(log_odds)

    @property
    def rank(self) -> int:
        return self.archetypes.shape[0]

    def design(self, features: np.ndarray):
        return structured.tensor_design(
            features, self.domain, self.weight_grid)

    def weights(self, design) -> np.ndarray:
        logits = np.zeros((design.shape[0], self.rank))
        for index in range(self.rank - 1):
            logits[:, index] = design @ self.log_odds[index]
        logits -= logits.max(axis=1, keepdims=True)
        weights = np.exp(logits)
        return weights / weights.sum(axis=1, keepdims=True)

    def shapes(self, features: np.ndarray) -> np.ndarray:
        if self.rank == 1:
            row = self.archetypes[0] / max(self.archetypes[0].sum(), 1e-300)
            return np.tile(row, (features.shape[0], 1))
        shapes = self.weights(self.design(features)) @ self.archetypes
        return shapes / np.maximum(
            shapes.sum(axis=1, keepdims=True), 1e-300)

    def normalize_rows(self) -> None:
        """Exact gauge fix: archetype rows to unit sum, the compensating
        constant folded into the log-odds fields."""
        sums = self.archetypes.sum(axis=1)
        floor = 1e-9 * max(float(sums.max()), 1e-300)
        safe = np.maximum(sums, floor)
        self.archetypes = self.archetypes / safe[:, None]
        if self.rank >= 2 and self.log_odds is not None:
            reference = float(np.log(safe[-1]))
            for index in range(self.rank - 1):
                self.log_odds[index] += float(np.log(safe[index])) - reference

    def copy(self) -> "ShapeModel":
        return ShapeModel(
            self.domain.copy(), self.weight_grid, self.archetypes.copy(),
            None if self.log_odds is None else self.log_odds.copy())


class JointFitter:
    """Alternating U / L minimization of the deployed error.

    `loss` selects the residual norm that stands in for total variation:
      l2        weighted least squares on the grouped residual
      cell_l1   sum_i c_i ||r_i||_2 by cell-level IRLS: robust to the few
                cells whose shape the mixture cannot represent at all, which
                is what the per-cell MEDIAN total variation rewards
      group_l1  sum_i c_i ||r_i||_1 by per-group IRLS
    """

    def __init__(self, data: ShapeData, model: ShapeModel,
                 weight_lambda: float, count_cap: float = 1.0e4,
                 count_weighting: bool = True, loss: str = "l2",
                 irls_floor: float = 0.02):
        if loss not in LOSSES:
            raise ValueError(f"unknown loss {loss!r}")
        self.data = data
        self.model = model
        self.rank = model.rank
        self.design = model.design(data.feature)
        self.field_size = self.design.shape[1]
        self.penalty = structured.second_difference_penalty(
            model.weight_grid, (weight_lambda,) * 3)
        counts = np.clip(data.counts, 1.0, count_cap)
        if not count_weighting:
            counts = np.ones_like(counts)
        self.base_weight = counts / counts.mean()
        self.cell_weight = self.base_weight.copy()
        self.loss = loss
        self.irls_floor = irls_floor
        self.group_weight = None
        self.matrix = data.matrix
        self.truth = data.truth
        self._rebuild_gram_pieces()

    def _rebuild_gram_pieces(self) -> None:
        scaled = (
            self.matrix if self.group_weight is None
            else self.matrix * self.group_weight[:, :, None])
        self.q_flat = np.einsum(
            "igk,igl->ikl", scaled, self.matrix).reshape(
                self.data.size, BASIS_SIZE * BASIS_SIZE)
        self.v = np.einsum("igk,ig->ik", scaled, self.truth)

    def reweight(self) -> None:
        """IRLS update of the residual weights from the current model."""
        if self.loss == "l2":
            return
        residual = self.residual()
        if self.loss == "cell_l1":
            norms = np.sqrt(np.sum(residual * residual, axis=1))
            self.cell_weight = self.base_weight / np.maximum(
                norms, self.irls_floor)
            self.cell_weight /= self.cell_weight.mean()
        else:
            self.group_weight = 1.0 / np.maximum(
                np.abs(residual),
                self.irls_floor / np.sqrt(residual.shape[1]))
            self.group_weight /= self.group_weight.mean()
            self._rebuild_gram_pieces()

    def model_shapes(self, weights: np.ndarray | None = None) -> np.ndarray:
        if weights is None:
            weights = self.model.weights(self.design)
        shapes = weights @ self.model.archetypes
        return shapes / np.maximum(
            shapes.sum(axis=1, keepdims=True), 1e-300)

    def residual(self, weights: np.ndarray | None = None) -> np.ndarray:
        predicted = np.einsum(
            "igk,ik->ig", self.matrix, self.model_shapes(weights))
        return predicted - self.truth

    def objective(self) -> float:
        """The surrogate being minimized (data term + field penalty)."""
        residual = self.residual()
        if self.loss == "group_l1":
            data_term = float(np.sum(
                self.cell_weight[:, None] * np.abs(residual)))
        elif self.loss == "cell_l1":
            data_term = float(np.sum(
                self.base_weight
                * np.sqrt(np.sum(residual * residual, axis=1))))
        else:
            data_term = float(np.sum(
                self.cell_weight * np.sum(residual * residual, axis=1)))
        penalty_term = 0.0
        if self.model.log_odds is not None:
            for row in self.model.log_odds:
                penalty_term += float(row @ (self.penalty @ row))
        return data_term + penalty_term

    def train_tv(self) -> np.ndarray:
        return deployed_tv(
            self.matrix, self.truth, self.model.shapes(self.data.feature))

    def update_archetypes(self, weights: np.ndarray,
                          unit_sum_rho: float = 1.0e4) -> None:
        rank = self.rank
        size = rank * BASIS_SIZE
        weighted = self.cell_weight[:, None] * weights
        pair = np.einsum("ir,is->irs", weighted, weights).reshape(
            self.data.size, rank * rank)
        gram = (pair.T @ self.q_flat).reshape(
            rank, rank, BASIS_SIZE, BASIS_SIZE)
        gram = np.transpose(gram, (0, 2, 1, 3)).reshape(size, size)
        rhs = (weighted.T @ self.v).reshape(size)
        # Unit-sum rows: rho * (1^T u_r - 1)^2 per archetype.
        scale = float(np.trace(gram)) / max(size, 1)
        rho = unit_sum_rho * max(scale, 1e-300)
        block = np.ones((BASIS_SIZE, BASIS_SIZE))
        for row in range(rank):
            lo = row * BASIS_SIZE
            gram[lo:lo + BASIS_SIZE, lo:lo + BASIS_SIZE] += rho * block
            rhs[lo:lo + BASIS_SIZE] += rho
        gram = 0.5 * (gram + gram.T)
        gram += 1.0e-9 * max(scale, 1e-300) * np.eye(size)
        factor = scipy.linalg.cholesky(gram, lower=False)
        target = scipy.linalg.solve_triangular(
            factor, rhs, lower=False, trans="T")
        solution, _ = scipy.optimize.nnls(factor, target, maxiter=50 * size)
        archetypes = solution.reshape(rank, BASIS_SIZE)
        keep = archetypes.sum(axis=1) > 1e-8
        if not np.all(keep):
            archetypes[~keep] = self.model.archetypes[~keep]
        self.model.archetypes = archetypes
        self.model.normalize_rows()

    def update_fields(self, iterations: int = 300) -> None:
        if self.rank < 2:
            return
        rank = self.rank
        size = self.field_size
        # Per-archetype grouped spectra, (n, G, R).
        projected = np.einsum(
            "igk,rk->igr", self.matrix, self.model.archetypes)
        penalty = self.penalty

        def cost(flat):
            fields = flat.reshape(rank - 1, size)
            logits = np.zeros((self.data.size, rank))
            for index in range(rank - 1):
                logits[:, index] = self.design @ fields[index]
            logits -= logits.max(axis=1, keepdims=True)
            weights = np.exp(logits)
            weights /= weights.sum(axis=1, keepdims=True)
            predicted = np.einsum("igr,ir->ig", projected, weights)
            residual = predicted - self.truth
            value = float(np.sum(
                self.cell_weight * np.sum(residual * residual, axis=1)))
            grad_w = 2.0 * self.cell_weight[:, None] * np.einsum(
                "ig,igr->ir", residual, projected)
            inner = np.sum(weights * grad_w, axis=1, keepdims=True)
            grad_z = weights * (grad_w - inner)
            gradient = np.empty((rank - 1, size))
            for index in range(rank - 1):
                gradient[index] = self.design.T @ grad_z[:, index]
                penalized = penalty @ fields[index]
                value += float(fields[index] @ penalized)
                gradient[index] += 2.0 * penalized
            return value, gradient.ravel()

        start = self.model.log_odds.ravel().copy()
        result = scipy.optimize.minimize(
            cost, start, jac=True, method="L-BFGS-B",
            options={"maxiter": iterations, "maxfun": 4 * iterations,
                     "ftol": 1e-12, "gtol": 1e-10})
        candidate = result.x.reshape(rank - 1, size)
        if (np.all(np.isfinite(candidate))
                and cost(result.x)[0] <= cost(start)[0]):
            self.model.log_odds = candidate

    def monitor(self) -> float:
        """Loss-independent progress signal: the weighted mean deployed
        total variation on the fitted rows.  Picks the best iterate, so an
        IRLS reweight (which moves the surrogate) can never leave the model
        worse than its initialization."""
        return float(np.average(self.train_tv(), weights=self.base_weight))

    def run(self, rounds: int = 8, tolerance: float = 1e-5) -> list[float]:
        history = [self.objective()]
        monitors = [self.monitor()]
        best = self.model.copy()
        best_monitor = monitors[0]
        for index in range(rounds):
            self.reweight()
            self.update_archetypes(self.model.weights(self.design))
            self.update_fields()
            history.append(self.objective())
            value = self.monitor()
            monitors.append(value)
            if value < best_monitor:
                best_monitor = value
                best = self.model.copy()
            if index > 0 and abs(monitors[-2] - value) <= tolerance * abs(
                    monitors[-2]):
                break
        self.model = best
        self.monitors = monitors
        return history


def pipeline_model(data: ShapeData, domain: np.ndarray, rank: int,
                   weight_grid, weight_lambda: float) -> ShapeModel:
    """The structured fitter's four-stage shape model on these rows: the
    joint fit's starting point."""
    stub = structured.ComponentData(
        identities=list(data.identity), features=data.feature,
        log_yield=data.log_yield, sigma_rel=data.sigma_rel,
        counts=data.counts, splits=data.split, shapes=data.shape)
    archetypes, fields, _ = structured.fit_shape_model(
        stub, np.ones(data.size, dtype=bool), domain, rank,
        tuple(weight_grid), weight_lambda)
    return ShapeModel(domain, weight_grid, archetypes, fields)


def load_settings(path: Path) -> dict:
    settings = json.loads(path.read_text())
    for name in SHAPE_COMPONENTS:
        spec = settings.get(name)
        if not isinstance(spec, dict):
            raise SystemExit(f"{path}: no settings for {name}")
        missing = [key for key in SETTING_KEYS if key not in spec]
        if missing:
            raise SystemExit(f"{path}: {name} lacks {', '.join(missing)}")
        grid = spec["weight_grid"]
        if (not 1 <= int(spec["rank"]) <= 16 or len(grid) != 3
                or not all(4 <= int(size) <= 32 for size in grid)
                or spec["loss"] not in LOSSES):
            raise SystemExit(
                f"{path}: {name} needs rank 1-16, a 3-axis weight grid "
                f"with sizes 4-32 and a loss in {LOSSES}")
    return settings


# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------

def command_extract(args) -> int:
    payload_path = args.fit_dir / PAYLOAD_NAME
    if not payload_path.exists():
        raise SystemExit(
            f"{payload_path} is missing: run the structured fit first")
    papf = structured.load_payload(
        payload_path).photon_absorption_path_fraction
    holdout_z = fit_holdout_z(args.fit_dir)
    shape_cache = structured.load_shape_cache(
        args.fit_dir / "cache" / f"shape_nnls_papf{papf:g}.pkl")
    basis_cache: dict = {}
    args.output.mkdir(parents=True, exist_ok=True)
    rows: dict[str, dict[str, int]] = {}
    # The off-grid FINAL TEST bank is never extracted: it is only scored.
    for label in ("bank", "offgrid_validation"):
        cache = args.fit_dir / "cache" / f"{CELL_CACHES[label]}.pkl"
        if label != "bank" and not cache.exists():
            log(f"no {label} bank in the structured fit; skipped")
            continue
        cells = fit_cells(args.fit_dir, label, holdout_z)
        log(f"{label}: {len(cells)} cells")
        rows[label] = {}
        for name in SHAPE_COMPONENTS:
            data = extract_component(
                cells, name, papf, shape_cache, basis_cache)
            np.savez_compressed(args.output / f"{label}__{name}.npz", **data)
            rows[label][name] = int(data["feature"].shape[0])
            log(f"  {name}: {rows[label][name]} rows")
    (args.output / "extract.json").write_text(json.dumps({
        "structured_payload_sha256": sha256(payload_path),
        "photon_absorption_path_fraction": papf,
        "holdout_z": sorted(holdout_z),
        "rows": rows}, indent=2) + "\n")
    return 0


# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------

def format_values(values) -> str:
    return " ".join(f"{float(value):.17g}" for value in values)


def shape_block_length(tokens: list[str]) -> int:
    rank = int(tokens[1])
    return 1 + rank + (rank if rank >= 2 else 0)


def shape_block_owner(tokens: list[str]) -> str | None:
    tag = tokens[0] if tokens else ""
    return next(
        (name for name in SHAPE_COMPONENTS if tag == f"{name}_shape_rank"),
        None)


def rewrite(source: Path, target: Path, models: dict) -> None:
    """Copy `source` with every shape block of a SHAPE_COMPONENT replaced
    by `models`; all other lines (header, domains, yield fields, the
    yield-only components) pass through unchanged."""
    lines = source.read_text().splitlines()
    output = []
    index = 0
    while index < len(lines):
        tokens = lines[index].split()
        name = shape_block_owner(tokens)
        if name is None:
            output.append(lines[index])
            index += 1
            continue
        index += shape_block_length(tokens)
        model = models[name]
        output.append(f"{name}_shape_rank {model.rank}")
        for row in model.archetypes:
            output.append(f"{name}_shape_basis " + format_values(row))
        if model.rank >= 2:
            output.append(
                f"{name}_weight_grid 3 " + " ".join(
                    str(size) for size in model.weight_grid))
            for row in model.log_odds:
                output.append(
                    f"{name}_weight_log_odds " + format_values(row))
    target.write_text("\n".join(output) + "\n")


def non_shape_lines(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    keep = []
    index = 0
    while index < len(lines):
        tokens = lines[index].split()
        if shape_block_owner(tokens) is None:
            keep.append(lines[index])
            index += 1
        else:
            index += shape_block_length(tokens)
    return keep


def command_build(args) -> int:
    settings = load_settings(args.settings)
    source = args.fit_dir / PAYLOAD_NAME
    extract_record = args.data / "extract.json"
    if not extract_record.exists():
        raise SystemExit(
            f"{extract_record} is missing: run the extract subcommand first")
    if (json.loads(extract_record.read_text())
            ["structured_payload_sha256"] != sha256(source)):
        raise SystemExit(
            f"the shape datasets in {args.data} were extracted from a "
            f"different structured fit than {source}: run extract again")
    reference = structured.load_payload(source)
    splits = sorted({value for value in args.fit_splits.split(",") if value})
    models, notes = {}, {}
    for name in SHAPE_COMPONENTS:
        spec = settings[name]
        bank = ShapeData.load(args.data / f"bank__{name}.npz")
        data = bank.subset(np.isin(bank.split, splits))
        if data.size == 0:
            raise SystemExit(
                f"{name}: no {'/'.join(splits)} rows to fit in {args.data}")
        start = time.time()
        model = pipeline_model(
            data, reference.components[name].domain, int(spec["rank"]),
            tuple(spec["weight_grid"]), float(spec["init_lambda"]))
        fitter = JointFitter(
            data, model, float(spec["weight_lambda"]),
            count_weighting=bool(spec["count_weighting"]),
            loss=spec["loss"])
        history = fitter.run(rounds=args.rounds)
        model = models[name] = fitter.model
        if not (np.all(model.archetypes >= 0.0) and np.allclose(
                model.archetypes.sum(axis=1), 1.0, rtol=0.0, atol=1e-12)):
            raise SystemExit(f"{name}: archetypes are not unit-sum rows")
        notes[name] = {
            "rows": int(data.size),
            **{key: spec[key] for key in SETTING_KEYS},
            "surrogate": [history[0], history[-1]],
            "deployed_tv_monitor": [fitter.monitors[0], min(fitter.monitors)],
            "rounds": len(history) - 1,
            "seconds": round(time.time() - start, 1)}
        log(f"{name}: {data.size} rows, R={model.rank}, deployed TV "
            f"{fitter.monitors[0]:.4f} -> {min(fitter.monitors):.4f}, "
            f"{time.time() - start:.0f} s")

    args.output.mkdir(parents=True, exist_ok=True)
    target = args.output / PAYLOAD_NAME
    rewrite(source, target, models)
    if non_shape_lines(source) != non_shape_lines(target):
        raise SystemExit("rewrite changed a non-shape payload line")
    parsed = structured.load_payload(target)
    rng = np.random.default_rng(7)
    for name in fit.MODEL_COMPONENTS:
        if name not in reference.components:
            continue   # a component the bank does not populate is "none"
        mine = parsed.components[name]
        base = reference.components[name]
        probe = rng.uniform(-0.2, 1.2, size=(64, 3)) * (
            mine.domain[:, 1] - mine.domain[:, 0]) + mine.domain[:, 0]
        if name in models:
            difference = np.abs(
                models[name].shapes(probe) - mine.predict_shape(probe))
            if float(difference.max()) > 1e-12:
                raise SystemExit(f"{name}: shape round trip failed")
        if not np.array_equal(
                mine.predict_yield(probe), base.predict_yield(probe)):
            raise SystemExit(f"{name}: yield field changed")
    log(f"payload written: {target} ({target.stat().st_size} bytes); "
        "non-shape lines identical to the structured payload, shape "
        "round trip < 1e-12")
    (args.output / "build_notes.json").write_text(json.dumps({
        "structured_payload_sha256": sha256(source),
        "payload_sha256": sha256(target),
        "fit_splits": splits, "rounds": args.rounds,
        "components": notes}, indent=2) + "\n")
    return 0


# ---------------------------------------------------------------------------
# score
# ---------------------------------------------------------------------------

METRICS = (
    "yield_relative_uncertainty", "yield_relative_error",
    "shape_total_variation", "weighted_rms_relative_error",
    "p95_relative_error")


def summarize_records(records: list[dict]) -> dict:
    answer = {}
    for name in fit.REPORT_COMPONENTS:
        selected = [record for record in records
                    if record["component"] == name]
        if selected:
            answer[name] = {
                metric: fit.finite_summary(
                    [abs(float(record[metric])) for record in selected])
                for metric in METRICS}
    return answer


def command_score(args) -> int:
    payload = structured.load_payload(args.payload)
    papf = payload.photon_absorption_path_fraction
    predict = structured.payload_predict_fn(payload)
    holdout_z = fit_holdout_z(args.fit_dir)
    wanted = [value for value in args.banks.split(",") if value]
    unknown = set(wanted) - {
        "test", "validation", "offgrid_validation", "offgrid_final_test"}
    if unknown:
        raise SystemExit(f"unknown banks: {', '.join(sorted(unknown))}")
    report = {
        "payload": args.payload.name,
        "payload_sha256": sha256(args.payload),
        "holdout_z": sorted(holdout_z),
        "banks": {}}
    basis_cache: dict = {}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    bank_cells = None
    for key in wanted:
        if key in ("test", "validation"):
            if bank_cells is None:
                bank_cells = fit_cells(args.fit_dir, "bank", holdout_z)
            cells = [cell for cell in bank_cells if cell.split == key]
        else:
            cells = fit_cells(args.fit_dir, key, holdout_z)
        if not cells:
            raise SystemExit(f"no {key} cells to score")
        log(f"scoring {key}: {len(cells)} cells")
        records, _ = fit.evaluate(cells, predict, basis_cache, papf)
        fit.write_metrics(args.output.parent / f"{key}_metrics.csv", records)
        report["banks"][key] = {
            "cells": len(cells), "components": summarize_records(records)}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    for key, block in report["banks"].items():
        log(f"\n== {key} ({block['cells']} cells): |yield error| and shape "
            "TV, median / p95")
        for name, metrics in block["components"].items():
            error = metrics["yield_relative_error"]
            tv = metrics["shape_total_variation"]
            if not error.get("count"):
                continue
            log(f"  {name:42s} n={error['count']:5d}  yield "
                f"{error['median']:7.2%} / {error['p95']:7.2%}   TV "
                f"{tv['median']:7.2%} / {tv['p95']:7.2%}")
    return 0


# ---------------------------------------------------------------------------
# parity
# ---------------------------------------------------------------------------

def choose_probes(cells: list, per_particle: int) -> list:
    """Evenly spaced cells, per particle, through the cells sorted by
    (energy, Z, areal density): deterministic and spanning the energy
    range."""
    chosen = []
    for positron in (False, True):
        group = sorted(
            (cell for cell in cells if cell.is_positron == positron),
            key=lambda cell: (cell.energy_keV, cell.z_value,
                              cell.areal_density, cell.identity))
        count = min(per_particle, len(group))
        if count == 0:
            continue
        picks = sorted({
            round(index * (len(group) - 1) / max(count - 1, 1))
            for index in range(count)})
        chosen.extend(group[index] for index in picks)
    return chosen


def cpp_spectra(evaluator: Path, payload: Path, cell, out: Path) -> dict:
    edges = np.asarray(cell.edges)
    command = [
        str(evaluator), "--payload", str(payload), "--output", str(out),
        "--energy-kev", repr(float(cell.energy_keV)),
        "--effective-z", repr(float(cell.z_value)),
        "--areal-density", repr(float(cell.areal_density)),
        "--bins", str(int(edges.size - 1)),
        "--max-energy-kev", repr(float(edges[-1]))]
    if cell.is_positron:
        command.append("--positron")
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        raise SystemExit(
            f"the C++ evaluator failed on {cell.identity} (exit "
            f"{completed.returncode}): {completed.stderr.strip()}\n  "
            + " ".join(command))
    answer: dict[str, list] = {}
    with out.open(newline="") as stream:
        for row in csv.DictReader(stream):
            answer.setdefault(row["component"], []).append(
                (float(row["bin_lo_keV"]), float(row["bin_hi_keV"]),
                 float(row["sum_w"])))
    return answer


def command_parity(args) -> int:
    if not os.access(args.evaluator, os.X_OK):
        raise SystemExit(f"{args.evaluator} is not an executable")
    payload = structured.load_payload(args.payload)
    papf = payload.photon_absorption_path_fraction
    predict = structured.payload_predict_fn(payload)
    holdout_z = fit_holdout_z(args.fit_dir)
    cells = [cell for cell in fit_cells(args.fit_dir, "bank", holdout_z)
             if cell.split == "test"]
    if not cells:
        raise SystemExit("the training bank has no whole-Z holdout cells")
    probes = choose_probes(cells, args.probes_per_particle)
    basis_cache: dict = {}
    rows = []
    with tempfile.TemporaryDirectory() as scratch:
        out = Path(scratch) / "cpp.csv"
        for cell in probes:
            blocks = cpp_spectra(args.evaluator, args.payload, cell, out)
            containment = None
            if cell.is_positron:
                containment = float(predict(
                    cell, "positron_containment", None, basis_cache,
                    papf)[0])
            for name in fit.MODEL_COMPONENTS:
                if (name == "positron_containment"
                        or not fit.component_applies(cell, name)):
                    continue  # containment is not a spectral output
                mine = predict(cell, name, containment, basis_cache, papf)
                entries = blocks.get(name, [])
                if name in fit.ATOM_COMPONENTS:
                    # Atoms (a fluorescence family's several lines) sum.
                    theirs = np.asarray([sum(row[2] for row in entries)])
                else:
                    # Match C++ bins to the cell's by bin centre: the two
                    # edge arrays can differ by an ulp at a boundary.
                    theirs = np.zeros(mine.size)
                    for lo, hi, value in entries:
                        index = int(np.searchsorted(
                            cell.edges, 0.5 * (lo + hi), "right") - 1)
                        theirs[min(max(index, 0), mine.size - 1)] += value
                scale = max(float(np.abs(mine).max()), 1e-300)
                difference = float(np.abs(theirs - mine).max()) / scale
                rows.append({"cell": cell.identity, "component": name,
                             "max_relative_difference": difference})
                log(f"  {cell.identity:34s} {name:42s} {difference:.2e}")
    worst = max(rows, key=lambda row: row["max_relative_difference"])
    args.output.write_text(json.dumps({
        "payload": args.payload.name,
        "payload_sha256": sha256(args.payload),
        "metric": "max over bins |C++ - Python| / max over bins |Python|",
        "probe_cells": len(probes),
        "worst_relative_difference": worst["max_relative_difference"],
        "worst_at": worst, "comparisons": rows}, indent=2) + "\n")
    log(f"worst relative difference over {len(probes)} probe cells: "
        f"{worst['max_relative_difference']:.3e} ({worst['cell']}, "
        f"{worst['component']})")
    return 0


# ---------------------------------------------------------------------------
# reselect (optional)
# ---------------------------------------------------------------------------

SELECT_MIN_SHAPE_COUNTS = 200.0
SCAN_LAMBDAS = (0.0, 1e-9, 1e-8, 1e-7, 1e-6)
SCAN_LOSSES = (("l2", True), ("l2", False), ("cell_l1", True),
               ("cell_l1", False))
SCAN_RANKS = (4, 6, 8, 12, 16)
SCAN_GRIDS = ((8, 6, 10), (10, 6, 12), (12, 8, 14))
PARSIMONY_SLACK = 0.02


def selection_score(model: ShapeModel, select) -> tuple[float, dict]:
    """Count-weighted selection objective over the selection sets, scored
    on rows with resolved truth (>= 200 counts) where there are any."""
    summaries, per_set = [], {}
    for label, data in select:
        if data.size == 0:
            continue
        mask = data.counts >= SELECT_MIN_SHAPE_COUNTS
        if not np.any(mask):
            mask = np.ones(data.size, dtype=bool)
        sub = data.subset(mask)
        summary = summarize(
            deployed_tv(sub.matrix, sub.truth, model.shapes(sub.feature)))
        per_set[label] = summary
        summaries.append((summary, summary["count"]))
    if not summaries:
        return math.inf, per_set
    return sum(selection_objective(summary) * count
               for summary, count in summaries) / sum(
                   count for _, count in summaries), per_set


def command_reselect(args) -> int:
    """Staged coordinate scan from the structured structure: weight-field
    smoothness, then residual norm and cell weighting, then rank, then
    weight grid (the rank and grid stages take the smallest candidate
    within 2% of the best)."""
    structure = json.loads(args.structure.read_text())
    components = [value for value in args.components.split(",") if value]
    settings, scan = {}, {}
    for name in components:
        spec = structure[name]
        bank = ShapeData.load(args.data / f"bank__{name}.npz")
        train = bank.subset(bank.split == "train")
        select = [("bank_validation", bank.subset(bank.split == "validation"))]
        offgrid_path = args.data / f"offgrid_validation__{name}.npz"
        if offgrid_path.exists():
            select.append(("offgrid_validation", ShapeData.load(offgrid_path)))
        if train.size == 0:
            raise SystemExit(f"{name}: no train rows in {args.data}")
        domain = structured.training_domain(train.feature)
        state = {"rank": int(spec["shape_rank"]),
                 "grid": tuple(spec["weight_grid"]), "lam": 0.0,
                 "loss": "l2", "count_weighting": True}
        tables = {}

        def trial(**overrides):
            trial_state = dict(state, **overrides)
            model = pipeline_model(
                train, domain, trial_state["rank"], trial_state["grid"],
                spec["weight_lambda"])
            fitter = JointFitter(
                train, model, trial_state["lam"],
                count_weighting=trial_state["count_weighting"],
                loss=trial_state["loss"])
            fitter.run(rounds=args.rounds)
            objective, per_set = selection_score(fitter.model, select)
            return {"objective": objective, "per_set": per_set}

        def smallest_within_slack(table, candidates, key):
            floor = min(table[key(c)]["objective"] for c in candidates)
            return next(
                c for c in candidates
                if table[key(c)]["objective"]
                <= floor + PARSIMONY_SLACK * max(floor, 0.01))

        pipeline = pipeline_model(
            train, domain, state["rank"], state["grid"],
            spec["weight_lambda"])
        start_objective, _ = selection_score(pipeline, select)
        table = {f"{lam:g}": trial(lam=lam) for lam in SCAN_LAMBDAS}
        state["lam"] = min(
            SCAN_LAMBDAS, key=lambda lam: table[f"{lam:g}"]["objective"])
        tables["weight_lambda"] = table
        table = {f"{loss}_{'count' if weighted else 'uniform'}": trial(
            loss=loss, count_weighting=weighted)
            for loss, weighted in SCAN_LOSSES}
        best = min(table, key=lambda key: table[key]["objective"])
        state["loss"], weighting = best.rsplit("_", 1)
        state["count_weighting"] = weighting == "count"
        tables["loss"] = table
        table = {str(rank): trial(rank=rank) for rank in SCAN_RANKS}
        state["rank"] = smallest_within_slack(table, SCAN_RANKS, str)
        tables["rank"] = table
        table = {str(list(grid)): trial(grid=grid) for grid in SCAN_GRIDS}
        state["grid"] = smallest_within_slack(
            table, sorted(SCAN_GRIDS, key=lambda grid: int(np.prod(grid))),
            lambda grid: str(list(grid)))
        tables["weight_grid"] = table
        final = trial()
        settings[name] = {
            "rank": state["rank"], "weight_grid": list(state["grid"]),
            "weight_lambda": state["lam"], "loss": state["loss"],
            "count_weighting": state["count_weighting"],
            "init_lambda": spec["weight_lambda"]}
        scan[name] = {"four_stage_objective": start_objective,
                      "selected_objective": final["objective"],
                      "selected": settings[name], "tables": tables}
        log(f"{name}: {settings[name]}  objective {start_objective:.4f} "
            f"(four-stage) -> {final['objective']:.4f}")
    args.output.write_text(json.dumps(settings, indent=2) + "\n")
    scan_path = args.output.with_name(args.output.stem + "_scan.json")
    scan_path.write_text(json.dumps(scan, indent=2) + "\n")
    log(f"settings written to {args.output}; scan tables in {scan_path}")
    return 0


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    fit_dir_help = "output directory of fit_charged_particle_source_structured.py"

    extract = commands.add_parser(
        "extract", help="per-component shape datasets from the structured "
        "fit's parsed-bank caches")
    extract.add_argument("--fit-dir", type=Path, required=True,
                         help=fit_dir_help)
    extract.add_argument("--output", type=Path, required=True,
                         help="directory for the <bank>__<component>.npz "
                         "datasets")

    build = commands.add_parser(
        "build", help="joint-fit the shape blocks and write the payload")
    build.add_argument("--fit-dir", type=Path, required=True,
                       help=fit_dir_help)
    build.add_argument("--data", type=Path, required=True,
                       help="the extract subcommand's --output")
    build.add_argument("--settings", type=Path, required=True,
                       help="per-component joint-fit settings (JSON)")
    build.add_argument("--output", type=Path, required=True,
                       help=f"directory for {PAYLOAD_NAME} and "
                       "build_notes.json")
    build.add_argument("--rounds", type=int, default=16,
                       help="alternation rounds per component")
    build.add_argument("--fit-splits", default="train,validation",
                       help="training-bank splits the shapes are fitted on")

    score = commands.add_parser(
        "score", help="score a payload on held-out banks")
    score.add_argument("--fit-dir", type=Path, required=True,
                       help=fit_dir_help)
    score.add_argument("--payload", type=Path, required=True)
    score.add_argument(
        "--banks", default="test,offgrid_validation,offgrid_final_test",
        help="comma list of test (the whole-Z holdout planes), validation, "
        "offgrid_validation, offgrid_final_test")
    score.add_argument("--output", type=Path, required=True,
                       help="report JSON; <bank>_metrics.csv files are "
                       "written next to it")

    parity = commands.add_parser(
        "parity", help="compare the C++ evaluator with the Python "
        "reference evaluator on whole-Z holdout cells")
    parity.add_argument("--fit-dir", type=Path, required=True,
                        help=fit_dir_help)
    parity.add_argument("--payload", type=Path, required=True)
    parity.add_argument("--evaluator", type=Path, required=True,
                        help="ceelo_charged_particle_model_evaluate")
    parity.add_argument("--probes-per-particle", type=int, default=6)
    parity.add_argument("--output", type=Path, required=True)

    reselect = commands.add_parser(
        "reselect", help="optional: re-select the joint-fit settings")
    reselect.add_argument("--data", type=Path, required=True,
                          help="the extract subcommand's --output")
    reselect.add_argument("--structure", type=Path, required=True,
                          help="the structured fit's structure JSON (start "
                          "rank, weight grid, initialization lambda)")
    reselect.add_argument("--components", default=",".join(SHAPE_COMPONENTS))
    reselect.add_argument("--rounds", type=int, default=12)
    reselect.add_argument("--output", type=Path, required=True,
                          help="settings JSON to write")

    args = parser.parse_args()
    return {
        "extract": command_extract, "build": command_build,
        "score": command_score, "parity": command_parity,
        "reselect": command_reselect,
    }[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
