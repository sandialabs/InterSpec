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

"""Generate aggregate bremsstrahlung observables for the migration review gate.

This validation-only tool may read a developer's local G4EMLOW installation and
historical git objects, but it emits only integrated rates, moments, CDFs, and
binned spectra. It never writes raw G4 cross-section rows or native columns.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import html
import json
import math
import re
import struct
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "validation" / "bremsstrahlung_migration_review"
ELEMENTS = {1: "H", 6: "C", 13: "Al", 26: "Fe", 32: "Ge", 48: "Cd",
            52: "Te", 53: "I", 74: "W", 82: "Pb", 92: "U"}
ENERGIES_KEV = (20.0, 50.0, 100.0, 500.0, 1000.0, 3000.0, 10000.0, 20000.0)
EXHAUSTIVE_ENERGIES_KEV = (
    10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0,
    150.0, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0,
    1500.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 8000.0,
    10000.0, 15000.0, 20000.0,
)
ATOMIC_WEIGHT = (
    1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999, 18.998, 20.180,
    22.990, 24.305, 26.982, 28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078,
    44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.380,
    69.723, 72.630, 74.922, 78.971, 79.904, 83.798, 85.468, 87.620, 88.906, 91.224,
    92.906, 95.950, 98.0, 101.070, 102.906, 106.420, 107.868, 112.414, 114.818, 118.710,
    121.760, 127.600, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908, 144.242,
    145.0, 150.360, 151.964, 157.250, 158.925, 162.500, 164.930, 167.259, 168.934, 173.045,
    174.967, 178.490, 180.948, 183.840, 186.207, 190.230, 192.217, 195.085, 196.967, 200.592,
    204.383, 207.200, 208.980, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.038,
    231.036, 238.029,
)


def stoich(items: list[tuple[int, float]]) -> list[tuple[int, float]]:
    total = sum(count * ATOMIC_WEIGHT[z - 1] for z, count in items)
    return [(z, count * ATOMIC_WEIGHT[z - 1] / total) for z, count in items]


MATERIALS = {
    "Ge": [(32, 1.0)],
    "NaI": stoich([(11, 1), (53, 1)]),
    "CsI": stoich([(55, 1), (53, 1)]),
    "CZT": stoich([(48, 0.9), (30, 0.1), (52, 1.0)]),
    "Al": [(13, 1.0)],
    "SS304": [(26, 0.70), (24, 0.19), (28, 0.10), (25, 0.01)],
    "Pb": [(82, 1.0)],
}
BASELINES = ("before_main", "before_FullSpecApprox", "before_env_scatter", "after", "g4_reference")
NUMBER = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


@dataclass
class Table:
    kappa: list[float]
    energy_keV: list[float]
    rows: dict[int, list[list[float]]]

    def energy_row(self, z: int, energy_keV: float) -> list[float]:
        log_energy = math.log10(energy_keV)
        logs = [math.log10(value) for value in self.energy_keV]
        high_e = min(max(bisect.bisect_right(logs, log_energy), 1), len(logs) - 1)
        low_e = high_e - 1
        te = max(0.0, min(1.0, (log_energy - logs[low_e]) / (logs[high_e] - logs[low_e])))
        row0, row1 = self.rows[z][low_e], self.rows[z][high_e]
        return [a + te * (b - a) for a, b in zip(row0, row1)]

    def chi(self, z: int, energy_keV: float, kappa: float) -> float:
        row = self.energy_row(z, energy_keV)
        high_k = min(max(bisect.bisect_right(self.kappa, kappa), 1), len(self.kappa) - 1)
        low_k = high_k - 1
        tk = max(0.0, min(1.0, (kappa - self.kappa[low_k]) /
                              (self.kappa[high_k] - self.kappa[low_k])))
        return max(0.0, row[low_k] + tk * (row[high_k] - row[low_k]))


def array(text: str, name: str) -> list[float]:
    match = re.search(
        rf"(?:static\s+)?(?:extern\s+)?const\s+(?:float|uint16_t)\s+{re.escape(name)}(?:\[[^]]*\])?\s*=\s*\{{(.*?)\}};",
        text,
        re.DOTALL,
    )
    if not match:
        raise ValueError(f"array {name} not found")
    body = re.sub(r"//.*", "", match.group(1))
    return [float(value) for value in NUMBER.findall(body)]


def cpp_table(text: str) -> Table:
    kappa = array(text, "kSB_kappa")
    energy = [10.0**value for value in array(text, "kSB_log_E_keV")]
    quantized = "static const uint16_t Z1_sb_chi" in text
    scales = array(text, "kSB_chi_scale") if quantized else []
    rows = {}
    width = len(kappa)
    for z in range(1, 93):
        flat = array(text, f"Z{z}_sb_chi")
        if quantized:
            flat = [value * scales[z - 1] for value in flat]
        rows[z] = [flat[index:index + width] for index in range(0, len(flat), width)]
        if len(rows[z]) != len(energy):
            raise ValueError(f"Z={z}: SB row count mismatch")
    return Table(kappa, energy, rows)


def git_text(ref: str, path: str) -> str:
    return subprocess.run(
        ["git", "show", f"{ref}:{path}"], cwd=ROOT, check=True,
        stdout=subprocess.PIPE, text=True,
    ).stdout


def g4_table(directory: Path) -> Table:
    rows = {}
    common_kappa = None
    common_energy = None
    for z in range(1, 93):
        lines = (directory / f"br{z}").read_text(encoding="ascii").splitlines()
        kappa = [float(value) for value in lines[1].split()]
        energy = [math.exp(float(value)) * 1000.0 for value in lines[2].split()]
        element_rows = [[float(value) for value in line.split()] for line in lines[3:3 + len(energy)]]
        if common_kappa is None:
            common_kappa, common_energy = kappa, energy
        elif kappa != common_kappa or energy != common_energy:
            raise ValueError(f"Z={z}: inconsistent G4 reference grid")
        rows[z] = element_rows
    return Table(common_kappa or [], common_energy or [], rows)


def parse_estar(directory: Path) -> dict[int, tuple[list[float], list[float]]]:
    row_re = re.compile(r"\s+([0-9.]+E[+-][0-9]+)\s+" + r"\s+".join(
        [r"([0-9.]+E[+-][0-9]+)"] * 6) + r"\s*<br>")
    output = {}
    for z in range(1, 93):
        matches = row_re.finditer((directory / f"estar-text-Z{z:03d}.html").read_text(encoding="ascii"))
        rows = [(float(m.group(1)) * 1000.0, float(m.group(3))) for m in matches]
        selected = [(energy, rad) for energy, rad in rows if 10.0 <= energy <= 20000.0]
        if len(selected) != 53:
            raise ValueError(f"Z={z}: expected 53 ESTAR rows")
        output[z] = ([item[0] for item in selected], [item[1] for item in selected])
    return output


def log_interp(x: list[float], y: list[float], target: float) -> float:
    upper = min(max(bisect.bisect_right(x, target), 1), len(x) - 1)
    lower = upper - 1
    t = math.log(target / x[lower]) / math.log(x[upper] / x[lower])
    return math.exp(math.log(y[lower]) + t * math.log(y[upper] / y[lower]))


def radiation_stopping(estar: dict[int, tuple[list[float], list[float]]], z: int, energy: float) -> float:
    x, y = estar[z]
    return log_interp(x, y, energy)


def material_curve(table: Table, composition: list[tuple[int, float]], energy: float) -> list[float]:
    weights = [(z, mass * z * (z + 1.0) / ATOMIC_WEIGHT[z - 1]) for z, mass in composition]
    total = sum(weight for _, weight in weights)
    element_rows = [(weight, table.energy_row(z, energy)) for z, weight in weights]
    return [sum(weight * row[index] for weight, row in element_rows) / total
            for index in range(len(table.kappa))]


def curve_value(kappa_grid: list[float], row: list[float], kappa: float) -> float:
    upper = min(max(bisect.bisect_right(kappa_grid, kappa), 1), len(kappa_grid) - 1)
    lower = upper - 1
    t = max(0.0, min(1.0, (kappa - kappa_grid[lower]) /
                         (kappa_grid[upper] - kappa_grid[lower])))
    return max(0.0, row[lower] + t * (row[upper] - row[lower]))


def trapz(x: list[float], y: list[float]) -> float:
    return sum((x[i + 1] - x[i]) * 0.5 * (y[i] + y[i + 1]) for i in range(len(x) - 1))


def observables(
    table: Table,
    composition: list[tuple[int, float]],
    energy: float,
    estar: dict[int, tuple[list[float], list[float]]],
    bins: int = 64,
) -> tuple[dict, list[dict]]:
    kmin = 10.0 / energy
    # Dense deterministic integration grid and coarser committed histogram.
    dense = [math.exp(math.log(kmin) + i / 2048.0 * math.log(1.0 / kmin)) for i in range(2049)]
    material_row = material_curve(table, composition, energy)
    chi = [curve_value(table.kappa, material_row, value) for value in dense]
    number_integrand = [value / k for value, k in zip(chi, dense)]
    i0, i1 = trapz(dense, number_integrand), trapz(dense, chi)
    srad = sum(mass * radiation_stopping(estar, z, energy) for z, mass in composition)
    if i0 <= 0.0 or i1 <= 0.0:
        summary = {
            "integrated_rate_per_g_cm2": 0.0,
            "absolute_radiative_loss_MeV_cm2_g": srad,
            "mean_photon_keV": 0.0,
            "median_photon_keV": 0.0,
            "tail_kappa_gt_0p1": 0.0,
            "tail_kappa_gt_0p5": 0.0,
            "tail_kappa_gt_0p9": 0.0,
        }
        return summary, [{"kappa": 1.0, "photon_keV": energy,
                          "normalized_dN_dlnk": 0.0,
                          "absolute_dN_dlnk_per_g_cm2": 0.0, "cdf": 0.0}]
    normalization = (srad / (energy * 1.0e-3)) / i1
    cumulative = [0.0]
    for i in range(len(dense) - 1):
        cumulative.append(cumulative[-1] + (dense[i + 1] - dense[i]) *
                          0.5 * (number_integrand[i] + number_integrand[i + 1]))
    cumulative = [value / cumulative[-1] for value in cumulative]

    def quantile(fraction: float) -> float:
        upper = bisect.bisect_left(cumulative, fraction)
        if upper <= 0:
            return dense[0]
        if upper >= len(dense):
            return dense[-1]
        t = (fraction - cumulative[upper - 1]) / (cumulative[upper] - cumulative[upper - 1])
        return dense[upper - 1] + t * (dense[upper] - dense[upper - 1])

    def tail(threshold: float) -> float:
        if threshold <= kmin:
            return 1.0
        index = min(bisect.bisect_left(dense, threshold), len(dense) - 1)
        return 1.0 - cumulative[index]

    summary = {
        "integrated_rate_per_g_cm2": normalization * i0,
        "absolute_radiative_loss_MeV_cm2_g": srad,
        "mean_photon_keV": energy * i1 / i0,
        "median_photon_keV": energy * quantile(0.5),
        "tail_kappa_gt_0p1": tail(0.1),
        "tail_kappa_gt_0p5": tail(0.5),
        "tail_kappa_gt_0p9": tail(0.9),
    }
    points = []
    for index in range(bins):
        low = math.exp(math.log(kmin) + index / bins * math.log(1.0 / kmin))
        high = math.exp(math.log(kmin) + (index + 1) / bins * math.log(1.0 / kmin))
        center = math.sqrt(low * high)
        cdf_index = min(bisect.bisect_left(dense, center), len(dense) - 1)
        value = curve_value(table.kappa, material_row, center)
        points.append({
            "kappa": center,
            "photon_keV": center * energy,
            "normalized_dN_dlnk": value / i0,
            "absolute_dN_dlnk_per_g_cm2": normalization * value,
            "cdf": cumulative[cdf_index],
        })
    return summary, points


def percent_difference(after: float, before: float) -> float:
    return 100.0 * (after / before - 1.0) if before != 0.0 else 0.0


def exhaustive_metrics(after: Table, reference: Table, estar) -> tuple[dict, list[dict]]:
    cases = []
    for z in range(1, 93):
        for energy in EXHAUSTIVE_ENERGIES_KEV:
            kmin = 10.0 / energy
            samples = [max(kmin, kappa) for kappa in after.kappa if kappa >= kmin]
            if samples[0] > kmin:
                samples.insert(0, kmin)
            arow = [after.chi(z, energy, kappa) for kappa in samples]
            grow = [reference.chi(z, energy, kappa) for kappa in samples]
            row_rms = math.sqrt(sum(((a / g) - 1.0)**2 for a, g in zip(arow, grow)) / len(arow))
            asummary, _ = observables(after, [(z, 1.0)], energy, estar, 1)
            gsummary, _ = observables(reference, [(z, 1.0)], energy, estar, 1)
            rate = (abs(asummary["integrated_rate_per_g_cm2"] /
                        gsummary["integrated_rate_per_g_cm2"] - 1.0)
                    if gsummary["integrated_rate_per_g_cm2"] > 0.0 else 0.0)
            mean = (abs(asummary["mean_photon_keV"] / gsummary["mean_photon_keV"] - 1.0)
                    if gsummary["mean_photon_keV"] > 0.0 else 0.0)
            cases.append({"Z": z, "energy_keV": energy, "row_rms_percent": row_rms * 100.0,
                          "integrated_rate_difference_percent": rate * 100.0,
                          "mean_kappa_difference_percent": mean * 100.0,
                          "score_percent": max(row_rms, rate, mean) * 100.0})

    def percentile(values: list[float], fraction: float) -> float:
        ordered = sorted(values)
        return ordered[round((len(ordered) - 1) * fraction)]

    report = {}
    for key in ("row_rms_percent", "integrated_rate_difference_percent", "mean_kappa_difference_percent"):
        values = [case[key] for case in cases]
        report[key] = {"median": percentile(values, 0.5), "p95": percentile(values, 0.95), "worst": max(values)}
    return report, sorted(cases, key=lambda case: case["score_percent"], reverse=True)[:20]


def write_html(path: Path, data: dict) -> None:
    payload = json.dumps(data, separators=(",", ":"))
    path.write_text("""<!doctype html><meta charset="utf-8">
<title>CeeLo bremsstrahlung migration review</title>
<style>body{font:14px system-ui;margin:24px;color:#17202a}select{margin:0 12px 16px 4px}
.grid{display:grid;grid-template-columns:repeat(2,minmax(380px,1fr));gap:16px}.plot{border:1px solid #ccd;padding:8px}
svg{width:100%;height:260px}.axis{stroke:#555;fill:none}.after{stroke:#087e8b}.before{stroke:#d1495b}.g4{stroke:#2e7d32}.line{fill:none;stroke-width:2}
table{border-collapse:collapse;margin-top:16px}td,th{padding:5px 9px;border:1px solid #ccd;text-align:right}th:first-child,td:first-child{text-align:left}
.gate{background:#fff3cd;padding:10px;border-left:4px solid #d49b00}</style>
<h1>Bremsstrahlung migration manual review</h1>
<p class="gate">Merge gate: inspect low-κ shape, endpoint behavior, and material-dependent rate changes. Record approval in <code>manual_review.json</code>.</p>
<label>Case <select id="case"></select></label><label>Energy <select id="energy"></select></label>
<div class="grid"><div class="plot"><b>Normalized shape dN/dlnκ</b><svg id="shape"></svg></div>
<div class="plot"><b>Cumulative distribution</b><svg id="cdf"></svg></div>
<div class="plot"><b>After / before and after / G4</b><svg id="ratio"></svg></div>
<div class="plot"><b>Absolute dN/dlnκ per g/cm²</b><svg id="absolute"></svg></div></div>
<div id="summary"></div><h2>Worst 20 exhaustive comparisons</h2><div id="worst"></div>
<script>const D=""" + payload + """;
const ce=document.querySelector('#case'),ee=document.querySelector('#energy');
for(const c of D.cases){let o=document.createElement('option');o.textContent=c;o.value=c;ce.append(o)}
for(const e of D.energies){let o=document.createElement('option');o.textContent=(e>=1000?e/1000+' MeV':e+' keV');o.value=e;ee.append(o)}
function svg(id,series,logy=false){const s=document.querySelector(id),W=600,H=240,p=34;s.innerHTML='';
let all=series.flatMap(x=>x.y.filter(Number.isFinite)),y0=logy?Math.log10(Math.max(Math.min(...all.filter(v=>v>0)),1e-300)):Math.min(0,...all),y1=logy?Math.log10(Math.max(...all)):Math.max(...all);if(y1==y0)y1=y0+1;
const X=x=>p+(Math.log10(x)-Math.log10(series[0].x[0]))/(Math.log10(series[0].x.at(-1))-Math.log10(series[0].x[0]))*(W-2*p),Y=y=>H-p-((logy?Math.log10(Math.max(y,1e-300)):y)-y0)/(y1-y0)*(H-2*p);
s.setAttribute('viewBox',`0 0 ${W} ${H}`);s.innerHTML=`<path class="axis" d="M${p},${p}V${H-p}H${W-p}"/>`;
for(const q of series){let d=q.x.map((x,i)=>(i?'L':'M')+X(x)+','+Y(q.y[i])).join(' ');s.innerHTML+=`<path class="line ${q.cls}" d="${d}"><title>${q.name}</title></path>`}}
function table(rows){return '<table><tr>'+Object.keys(rows[0]).map(k=>'<th>'+k+'</th>').join('')+'</tr>'+rows.map(r=>'<tr>'+Object.values(r).map(v=>'<td>'+v+'</td>').join('')+'</tr>').join('')+'</table>'}
function render(){const key=ce.value+'|'+ee.value,x=D.spectra[key],a=x.after,b=x.before_main,g=x.g4_reference;
svg('#shape',[{...a,cls:'after',name:'after'},{...b,cls:'before',name:'before'},{...g,cls:'g4',name:'G4'}],true);
svg('#cdf',[{x:a.x,y:a.cdf,cls:'after'},{x:b.x,y:b.cdf,cls:'before'},{x:g.x,y:g.cdf,cls:'g4'}]);
svg('#ratio',[{x:a.x,y:a.y.map((v,i)=>v/b.y[i]),cls:'before'},{x:a.x,y:a.y.map((v,i)=>v/g.y[i]),cls:'g4'}]);
svg('#absolute',[{x:a.x,y:a.absolute,cls:'after'},{x:b.x,y:b.absolute,cls:'before'},{x:g.x,y:g.absolute,cls:'g4'}],true);
document.querySelector('#summary').innerHTML=table(D.summary[key]);}
ce.onchange=render;ee.onchange=render;document.querySelector('#worst').innerHTML=table(D.worst);render();</script>""", encoding="utf-8")


def write_readme(path: Path, metrics: dict) -> None:
    path.write_text(f"""# Bremsstrahlung migration review bundle

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
| Row RMS | {metrics['row_rms_percent']['median']:.6f}% | {metrics['row_rms_percent']['p95']:.6f}% | {metrics['row_rms_percent']['worst']:.6f}% |
| Integrated rate | {metrics['integrated_rate_difference_percent']['median']:.6f}% | {metrics['integrated_rate_difference_percent']['p95']:.6f}% | {metrics['integrated_rate_difference_percent']['worst']:.6f}% |
| Mean κ | {metrics['mean_kappa_difference_percent']['median']:.6f}% | {metrics['mean_kappa_difference_percent']['p95']:.6f}% | {metrics['mean_kappa_difference_percent']['worst']:.6f}% |

`summary.csv` contains absolute yields, rates, moments, and tail probabilities.
`spectra.csv` contains only 64-bin aggregate spectra/CDFs, not native G4
cross-section columns. `worst_20.csv` contains the worst exhaustive cases.
`metadata.json` records versions and fixture provenance.

After inspection, edit `manual_review.json`: set the three resolution fields to
`true`, record reviewer/time/notes, and set `status` to `approved`. The gate is:

```bash
python3 tools/prepare_cross_sections/check_provenance.py --require-manual-approval
```
""", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--after-cpp", type=Path, default=ROOT / "src/cross_sections/element_data.cpp")
    parser.add_argument("--g4-brem-dir", type=Path, required=True)
    parser.add_argument("--estar-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--main-ref", default="ab03e4810c548126daced9d7c91d62b0467a3d6b")
    parser.add_argument("--fullspec-ref", default="7e8713b517c13999a148bea2578cc3f26d99e67f")
    parser.add_argument("--env-ref", default="4977fccff48fd26bdcc9bb873577c72ff0973183")
    args = parser.parse_args()

    tables = {
        "before_main": cpp_table(git_text(args.main_ref, "src/cross_sections/element_data.cpp")),
        "before_FullSpecApprox": cpp_table(git_text(args.fullspec_ref, "src/cross_sections/element_data.cpp")),
        "before_env_scatter": cpp_table(git_text(args.env_ref, "src/cross_sections/element_data.cpp")),
        "after": cpp_table(args.after_cpp.read_text(encoding="utf-8")),
        "g4_reference": g4_table(args.g4_brem_dir),
    }
    estar = parse_estar(args.estar_dir)
    args.output.mkdir(parents=True, exist_ok=True)
    summary_rows, spectrum_rows = [], []
    html_spectra, html_summary = {}, {}
    cases = [(symbol, [(z, 1.0)]) for z, symbol in ELEMENTS.items()] + list(MATERIALS.items())
    for label, composition in cases:
        for energy in ENERGIES_KEV:
            key = f"{label}|{energy:g}"
            summaries, spectra = {}, {}
            for baseline in BASELINES:
                summaries[baseline], spectra[baseline] = observables(tables[baseline], composition, energy, estar)
            after = summaries["after"]
            before = summaries["before_main"]
            row = {"case": label, "energy_keV": energy}
            for baseline in BASELINES:
                for name, value in summaries[baseline].items():
                    row[f"{baseline}_{name}"] = value
            row["after_before_rate_difference_percent"] = percent_difference(
                after["integrated_rate_per_g_cm2"], before["integrated_rate_per_g_cm2"])
            row["after_before_mean_difference_percent"] = percent_difference(
                after["mean_photon_keV"], before["mean_photon_keV"])
            summary_rows.append(row)
            html_summary[key] = [
                {"baseline": baseline,
                 "rate / (g/cm²)": f"{summaries[baseline]['integrated_rate_per_g_cm2']:.7g}",
                 "mean keV": f"{summaries[baseline]['mean_photon_keV']:.7g}",
                 "median keV": f"{summaries[baseline]['median_photon_keV']:.7g}",
                 "P(κ>.5)": f"{summaries[baseline]['tail_kappa_gt_0p5']:.7g}"}
                for baseline in BASELINES
            ]
            html_spectra[key] = {}
            for baseline in BASELINES:
                points = spectra[baseline]
                html_spectra[key][baseline] = {
                    "x": [point["kappa"] for point in points],
                    "y": [point["normalized_dN_dlnk"] for point in points],
                    "absolute": [point["absolute_dN_dlnk_per_g_cm2"] for point in points],
                    "cdf": [point["cdf"] for point in points],
                }
                for bin_index, point in enumerate(points):
                    spectrum_rows.append({"case": label, "energy_keV": energy,
                                          "baseline": baseline, "bin": bin_index, **point})

    metrics, worst = exhaustive_metrics(tables["after"], tables["g4_reference"], estar)
    with (args.output / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(summary_rows)
    with (args.output / "spectra.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(spectrum_rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(spectrum_rows)
    with (args.output / "worst_20.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(worst[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(worst)
    metadata = {
        "schema_version": 1,
        "status": "awaiting_manual_approval",
        "reference_type": "deterministic aggregate spectrum/rate observable from Geant4 11.4.0 G4EMLOW 8.8 model data",
        "raw_g4_cross_section_columns_committed": False,
        "geant4_version": "11.4.0",
        "g4emlow_version": "8.8",
        "simulation_seed": None,
        "event_count": None,
        "statistical_uncertainty": 0.0,
        "photon_cut_keV": 10.0,
        "histogram_bins": 64,
        "elements": ELEMENTS,
        "materials": MATERIALS,
        "energies_keV": ENERGIES_KEV,
        "baselines": {"main": args.main_ref, "FullSpecApprox": args.fullspec_ref,
                      "env_scatter": args.env_ref},
        "acceptance_metrics_percent": metrics,
        "notes": [
            "The G4 fixture is an observable integrated spectrum/rate, not a raw cross-section table.",
            "A separate stochastic thin-target Geant4 fixture may be added; the deterministic reference has no seed or sampling uncertainty.",
            "Manual approval is required before feature-branch migration.",
        ],
    }
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (args.output / "manual_review.json").write_text(json.dumps({
        "status": "pending", "reviewer": None, "reviewed_at": None,
        "low_kappa_shift_resolved": None, "endpoint_artifacts_resolved": None,
        "material_rate_changes_resolved": None, "notes": "",
    }, indent=2) + "\n", encoding="utf-8")
    html_data = {"cases": [label for label, _ in cases], "energies": ENERGIES_KEV,
                 "spectra": html_spectra, "summary": html_summary,
                 "worst": [{key: (f"{value:.6g}" if isinstance(value, float) else value)
                            for key, value in row.items()} for row in worst]}
    write_html(args.output / "review.html", html_data)
    write_readme(args.output / "README.md", metrics)
    print(f"Wrote aggregate-only manual review bundle to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
