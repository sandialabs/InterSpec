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

# Checks the machine before an unattended bake (all.sh runs it first): the work
# area's file system and free space, memory, Python, power and priority.  Hard
# problems exit 1; the rest are printed as warnings.
#
#   tools/full_spectrum/bake/preflight.sh
#
# FAMILIES (as for all.sh) sets the disk space it asks for; ENV_ROUNDS > 0 adds
# nothing here, because env_scatter.sh checks each round's space itself.
. "$(dirname "$0")/common.sh"
FAMILIES=${FAMILIES:-"detector source_escape charged_particle env_scatter"}
DETECTORS=${DETECTORS:-"hpge_coax"}

problems=0
problem() { printf 'PROBLEM: %s\n' "$*"; problems=$(( problems + 1 )); }
warning() { printf 'warning: %s\n' "$*"; }
ok() { printf 'ok:      %s\n' "$*"; }

mkdir -p "$WORK"
mount_point=$(df -P "$WORK" | awk 'NR == 2 { print $6 }')
fs_type=$(mount | awk -v m="$mount_point" '$3 == m { sub(/^\(/, "", $4); sub(/,$/, "", $4); print $4 }')
case $fs_type in
    apfs|hfs|ext4|xfs|btrfs|zfs) ok "WORK=$WORK is on $fs_type ($mount_point)" ;;
    "") warning "could not tell the file system of $WORK ($mount_point)" ;;
    *) problem "WORK=$WORK is on $fs_type: macOS writes a ._ AppleDouble file beside every file that carries an extended attribute on non-native volumes, and the bake tools would read them as data; use an APFS volume (copy results to other disks with cp -X)" ;;
esac
if find "$WORK" -name '._*' -print -quit 2>/dev/null | grep -q .; then
    problem "AppleDouble ._ files already exist under $WORK; find out where they came from before baking"
fi

# remaining DIR GIB: what a work area that reaches GIB still needs.
remaining() {
    local used=0
    [ -d "$1" ] && used=$(du -sk "$1" | awk '{ printf "%d", $1 / 1048576 }')
    echo $(( $2 > used ? $2 - used : 0 ))
}
# GiB each family's work area reaches (banks, fits, archive), with room to
# spare over what a full bake used.
need=20   # margin
for family in $FAMILIES; do
    case $family in
        detector)
            for detector in $DETECTORS; do
                need=$(( need + $(remaining "$WORK/detector_$detector" 1) ))
            done ;;
        source_escape) need=$(( need + $(remaining "$WORK/source_escape" 2) )) ;;
        charged_particle) need=$(( need + $(remaining "$WORK/charged_particle" 2) )) ;;
        # Baseline archive, folded adaptive rounds and the spool.
        env_scatter) need=$(( need + $(remaining "$WORK/env_field_soil" 35) )) ;;
    esac
done
[ "$SMOKE" = 1 ] && need=5   # smoke banks are megabytes
free=$(df -Pk "$WORK" | awk 'NR == 2 { printf "%d", $4 / 1048576 }')
if [ "$free" -lt "$need" ]; then
    problem "$free GiB free on $mount_point; the selected families need about $need GiB"
else
    ok "$free GiB free on $mount_point (the selected families need about $need GiB)"
fi

memory=$(( $(sysctl -n hw.memsize 2>/dev/null || awk '/MemTotal/ { print $2 * 1024 }' /proc/meminfo) / 1073741824 ))
case " $FAMILIES " in
    *" env_scatter "*)
        if [ "$memory" -lt 16 ]; then
            problem "$memory GiB of memory; environment scatter needs 16 GiB"
        elif [ "$memory" -lt 24 ]; then
            warning "$memory GiB of memory: environment scatter peaks near 16 GiB, so close other applications"
        else
            ok "$memory GiB of memory"
        fi ;;
    *) ok "$memory GiB of memory" ;;
esac

if "$PY" -c 'import sys, numpy, scipy, h5py; sys.exit(sys.version_info < (3, 10))' 2>/dev/null; then
    ok "Python $("$PY" -c 'import sys; print(sys.version.split()[0])') at $PY has numpy, scipy and h5py"
else
    problem "$PY is not Python >= 3.10 with numpy, scipy and h5py (pip install -r tools/full_spectrum/requirements.txt)"
fi

niceness=$(ps -o ni= -p $$ | tr -d ' ')
if [ "$niceness" = 20 ]; then
    ok "running at nice 20"
else
    warning "running at nice $niceness: start the bake with 'nice -n 20' so it yields to other work"
fi
if command -v pmset > /dev/null; then
    if pmset -g batt | grep -q "AC Power"; then
        ok "on AC power"
    else
        warning "on battery: connect power, or the machine may sleep or throttle"
    fi
    # `caffeinate -is UTILITY` execs UTILITY in its own process and holds the
    # assertion in a child, so look for a caffeinate whose parent is this
    # script or one of its ancestors.
    pid=$$ ancestors=" "
    while [ "${pid:-1}" -gt 1 ]; do
        ancestors="$ancestors$pid "
        pid=$(ps -o ppid= -p "$pid" | tr -d ' ')
    done
    caffeinated=0
    for parent in $(ps -axo ppid=,comm= | awk '$2 ~ /caffeinate$/ { print $1 }'); do
        case $ancestors in *" $parent "*) caffeinated=1 ;; esac
    done
    if [ "$caffeinated" = 1 ]; then
        ok "running under caffeinate"
    else
        warning "not running under caffeinate: start the bake with 'caffeinate -is' so the machine does not sleep"
    fi
fi
if command -v tmutil > /dev/null && tmutil status 2>/dev/null | grep -q 'Running = 1'; then
    warning "a Time Machine backup is running; it competes for the disk"
fi

if [ "$problems" -gt 0 ]; then
    echo "PREFLIGHT FAILED: $problems problem(s)"
    exit 1
fi
echo "PREFLIGHT OK"
