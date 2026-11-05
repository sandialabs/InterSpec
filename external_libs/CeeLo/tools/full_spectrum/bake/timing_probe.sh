#!/usr/bin/env bash
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

# Times a few representative pieces of bake work on this machine (about 5
# minutes), so the reference times in README.md can be rescaled: divide each
# time here by the reference time for the machine's speed factor.  Run it on
# an otherwise idle machine, after build.sh.
#
#   tools/full_spectrum/bake/timing_probe.sh
. "$(dirname "$0")/common.sh"
precision_tiers   # nothing to scale
family timing_probe
REDO=all
need_bin ceelo_detector_response_campaign ceelo_source_escape_campaign \
    ceelo_source_electron_campaign env_kernel_bake

host() {
    uname -a
    echo "cores: $(getconf _NPROCESSORS_ONLN 2>/dev/null)"
    grep -m1 'model name' /proc/cpuinfo 2>/dev/null \
        || sysctl -n machdep.cpu.brand_string 2>/dev/null || true
}

# timed LABEL OUTPUT COMMAND...: run COMMAND with its output in OUTPUT and
# print the wall time.
timed() {
    local label=$1 output=$2 start end
    shift 2
    start=$("$PY" -c 'import time; print(time.time())')
    "$@" > "$output" 2>&1 || die "$label failed; see $output"
    end=$("$PY" -c 'import time; print(time.time())')
    "$PY" -c "print(f'{\"$label\":58s} {$end - $start:7.1f} s wall')"
}

probe() {
    timed "detector cell 661.657 keV, 500k events, 1 thread" "$FW/det_662.csv" \
        "$BIN/ceelo_detector_response_campaign" --detector hpge_coax \
        --energy-kev 661.657 --cos-theta 1 --bins 256 \
        --min-events 500000 --max-events 500000 --seed 12345 \
        --detector-brems --detector-fluorescence --detector-spectrum-fidelity
    timed "detector cell 2614.511 keV, 3M events, 1 thread" "$FW/det_2614.csv" \
        "$BIN/ceelo_detector_response_campaign" --detector hpge_coax \
        --energy-kev 2614.511 --cos-theta 1 --bins 256 \
        --min-events 3000000 --max-events 3000000 --seed 12346 \
        --detector-brems --detector-fluorescence --detector-spectrum-fidelity
    timed "source-escape cell 3000 keV Z=26 tau=4, 1M events, 1 thread" "$FW/src_3000.csv" \
        "$BIN/ceelo_source_escape_campaign" --energy-kev 3000 --z 26 --tau 4 \
        --events 1000000 --bins 256 --seed 12347 --source-compton-doppler
    timed "electron cell 5000 keV Z=82 1 g/cm2, 60k events, 1 thread" "$FW/electron.log" \
        "$BIN/ceelo_source_electron_campaign" --energy-kev 5000 --z 82 \
        --areal-density 1 --events 60000 --seed 12348 --bins 512 \
        --max-energy-kev 5000 --min-photon-kev 10 --output "$FW/electron.csv"
    local energy
    for energy in 2614 59.5; do
        rm -rf "$FW/env_$energy"
        mkdir -p "$FW/env_$energy"
        timed "env node $energy keV at 5 m, 4 x 262144 histories, all cores" \
            "$FW/env_$energy.log" \
            "$BIN/env_kernel_bake" --mode field --material soil \
            --energies "$energy" --distances 500 --histories 262144 \
            --registered 384 --e-floor 10 --mu-bins 64 --brems-floor 150 \
            --max-vertices 50 --scoring-radius 4 --role train \
            --h-src 100 --h-det 100 --block-id 1 --block-count 4 \
            --fit-block "$FW/env_$energy/block_{block_id}.ceefit"
        rm -rf "$FW/env_$energy"
    done
}

step host host
step probe probe
stamp "TIMING_PROBE_DONE: compare with the reference table in tools/full_spectrum/README.md"
