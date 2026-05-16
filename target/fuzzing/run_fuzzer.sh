#!/usr/bin/env bash
# Helper script to run a fuzzer without polluting the committed seed corpus.
#
# Usage: ./run_fuzzer.sh <fuzzer-name> [extra libFuzzer args...]
#
# The seed corpus at corpus/<fuzzer-name>/ is read-only.  libFuzzer writes new
# coverage-expanding inputs to a temp directory created under $TMPDIR (or /tmp),
# which is removed on exit.  Crash artifacts are written to the temp dir and
# copied to ./crashes/ on a non-zero exit so they survive cleanup.
#
# Examples:
#   ./run_fuzzer.sh fuzz_color_theme -max_total_time=300
#   ./run_fuzzer.sh fuzz_app_url -max_total_time=600 -rss_limit_mb=4096

set -u

if [ $# -lt 1 ]; then
  echo "usage: $0 <fuzzer-name> [extra libFuzzer args...]" >&2
  exit 1
fi

FUZZER_NAME="$1"
shift

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
SEED_DIR="${SCRIPT_DIR}/corpus/${FUZZER_NAME#fuzz_}"
BIN="${SCRIPT_DIR}/build/${FUZZER_NAME}"

if [ ! -x "${BIN}" ]; then
  echo "error: ${BIN} not found - did you cmake --build target/fuzzing/build ?" >&2
  exit 1
fi
if [ ! -d "${SEED_DIR}" ]; then
  echo "error: seed dir ${SEED_DIR} not found" >&2
  exit 1
fi

WORK_DIR="$( mktemp -d -t "interspec_${FUZZER_NAME}_XXXXXX" )"
CRASH_OUT="${SCRIPT_DIR}/crashes/${FUZZER_NAME}"
mkdir -p "${CRASH_OUT}"

cleanup() {
  # Save any crash/leak/timeout artifacts before nuking the work dir.
  shopt -s nullglob
  local saved=0
  for f in "${WORK_DIR}"/crash-* "${WORK_DIR}"/leak-* "${WORK_DIR}"/timeout-*; do
    cp -f "$f" "${CRASH_OUT}/" && saved=$((saved+1))
  done
  if [ "${saved}" -gt 0 ]; then
    echo
    echo "=== Saved ${saved} artifact(s) to ${CRASH_OUT}/"
  fi
  rm -rf "${WORK_DIR}"
}
trap cleanup EXIT

# Default ASAN_OPTIONS unless caller already set them.  detect_container_overflow=0
# silences a libc++/libFuzzer false-positive in libFuzzer's corpus loader.
# detect_leaks=0 ignores per-iteration leaks (relevant for fuzz_app_url, where
# every Wt widget is "leaked" by design across iterations).
: "${ASAN_OPTIONS:=detect_container_overflow=0:detect_leaks=0:abort_on_error=1}"
export ASAN_OPTIONS

echo "=== ${FUZZER_NAME}"
echo "    seed: ${SEED_DIR}"
echo "    work: ${WORK_DIR}"
echo

# libFuzzer reads from BOTH dirs; only writes to the FIRST (work dir).
"${BIN}" -artifact_prefix="${WORK_DIR}/" -print_final_stats=1 "$@" "${WORK_DIR}" "${SEED_DIR}"
