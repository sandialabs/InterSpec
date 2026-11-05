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

"""Convert a CeeLo detector geometry from its XML form to the plain-text form.

Reads the <Detector> element that CeeLo's XML codec writes -- a standalone
<CeeLoGeometry> file, or an InterSpec DRF (<DetectorPeakResponse>), which
carries it under <CeeLoGeometry> and <CeeLoResponse> -- and writes
GeometryDescriptor's "CEELO_GEOMETRY 1" text, which the full-spectrum bake
tools read without rapidxml.  Attribute names and defaults follow
read_detector_node() in src/io/DetectorResponse.cpp.

Usage:
    drf_geometry.py INPUT.xml [--output FILE]
"""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def find_detector(root: ET.Element) -> ET.Element:
    """The <Detector> of a CeeLoGeometry block, else of a CeeLoResponse."""
    for path in ("CeeLoGeometry/Detector", ".//CeeLoGeometry/Detector",
                 ".//CeeLoResponse/Detector"):
        found = root.find(path)
        if found is not None:
            return found
    if root.tag == "CeeLoGeometry" and root.find("Detector") is not None:
        return root.find("Detector")
    if root.tag == "Detector":
        return root
    raise SystemExit("no <CeeLoGeometry>/<Detector> or <CeeLoResponse>/<Detector> "
                     "element in the input")


def number(element: ET.Element, name: str, default: float) -> float:
    value = element.get(name)
    return default if value is None else float(value)


def flag(element: ET.Element, name: str) -> bool:
    value = element.get(name)
    return bool(value) and value[0] in "1tTyY"


def text(value: float) -> str:
    return repr(float(value))   # shortest exact round trip


def to_geometry_text(detector: ET.Element) -> str:
    lines = ["CEELO_GEOMETRY 1"]
    lines.append("shape " + ("box" if detector.get("shape") == "box" else "cylinder"))
    dimensions = (detector.findtext("Dimensions") or "").split()
    if not dimensions:
        raise SystemExit("<Detector> has no <Dimensions>")
    lines.append("dimensions " + " ".join(text(float(v)) for v in dimensions))
    lines.append("bullet_radius " + text(number(detector, "bulletRadius", 0.0)))
    lines.append(f"crystal_material {int(number(detector, 'crystalMaterial', -1))}")
    bore = detector.find("Bore")
    if bore is not None:
        lines.append(f"bore {text(number(bore, 'radius', 0.0))} "
                     f"{text(number(bore, 'depth', 0.0))} "
                     f"{1 if flag(bore, 'roundedTip') else 0}")
    dead = detector.find("DeadLayer")
    if dead is not None:
        lines.append("dead_layer " + " ".join(
            text(number(dead, side, 0.0)) for side in ("front", "side", "back")))
    for layer in detector.findall("Layer"):
        lines.append(f"layer {int(number(layer, 'material', -1))} "
                     + " ".join(text(number(layer, name, 0.0)) for name in (
                         "frontThickness", "sideThickness", "zStart", "zEnd")))
    collimator = detector.find("Collimator")
    if collimator is not None:
        lines.append(f"collimator {int(number(collimator, 'material', -1))} "
                     + " ".join(text(number(collimator, name, 0.0)) for name in (
                         "sideThickness", "zStart", "zEnd")))
    lines.append("symmetry " + ("quadrant" if detector.get("symmetry") == "quadrant"
                                else "axial"))
    lines.append("reference_point " + (
        "endcap_front" if detector.get("referencePoint") == "endcap_front"
        else "crystal_face"))
    materials = detector.find("Materials")
    entries = [] if materials is None else materials.findall("Material")
    entries.sort(key=lambda m: int(number(m, "index", -1)))
    for expected, material in enumerate(entries):
        index = int(number(material, "index", -1))
        if index != expected:
            raise SystemExit(f"material indices must be 0, 1, ...: found {index} "
                             f"where {expected} was expected")
        elements = material.findall("El")
        fields = [str(index), text(number(material, "density", 0.0)), str(len(elements))]
        for element in elements:
            fields += [str(int(number(element, "Z", 0))), text(number(element, "frac", 0.0))]
        name = material.get("name") or f"material{index}"
        lines.append("material " + " ".join(fields) + " " + name)
    lines.append("end")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("input", type=Path, help="DRF or CeeLoGeometry XML file")
    parser.add_argument("--output", type=Path,
                        help="where to write the text (default: standard output)")
    args = parser.parse_args()
    detector = find_detector(ET.parse(args.input).getroot())
    result = to_geometry_text(detector)
    if args.output is None:
        sys.stdout.write(result)
    else:
        args.output.write_text(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
