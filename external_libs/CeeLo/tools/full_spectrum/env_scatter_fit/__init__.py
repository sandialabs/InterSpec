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

"""Fit-first environmental-scatter campaign utilities.

The design, allocation, block metadata, and conservative rebinning modules use
only the Python standard library.  Numerical fitting requires NumPy; archive
conversion additionally requires h5py.  Those imports are deliberately lazy so
that worker and campaign hosts do not inherit either dependency.
"""

from .axes import (
    AxisDefinition,
    annihilation_registered_axis,
    brems_registered_axis,
    conservative_rebin_2d,
    primary_registered_axis,
    registered_edges_from_physical,
    uniform_mu_edges,
)
from .calibration import (
    CALIBRATION_SCHEMA,
    analyze_uncertainty_calibration,
    calibrate_uncertainty_archive,
)
from .design import (
    DesignBounds,
    DesignSelection,
    PhysicalNode,
    generate_candidate_pool,
    geometry_transform,
    select_disjoint_maximin,
    transformed_domain,
)

__all__ = [
    "AxisDefinition",
    "CALIBRATION_SCHEMA",
    "DesignBounds",
    "DesignSelection",
    "PhysicalNode",
    "annihilation_registered_axis",
    "analyze_uncertainty_calibration",
    "brems_registered_axis",
    "calibrate_uncertainty_archive",
    "conservative_rebin_2d",
    "generate_candidate_pool",
    "geometry_transform",
    "primary_registered_axis",
    "registered_edges_from_physical",
    "select_disjoint_maximin",
    "transformed_domain",
    "uniform_mu_edges",
]

__version__ = "0.1.0"
