#!/bin/bash
# Creates the interspec-assets.zip that gets bundled into the APK.
# This zip includes InterSpec_resources/, data/, example_spectra/, and Wt resources/.
#
# Wt resources location is resolved from (in order):
#   1. WT_RESOURCES_DIR environment variable
#   2. Common local build directories (build_xcode, build_vscode, build)
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ASSETS_ZIP="${SCRIPT_DIR}/app/src/main/assets/interspec-assets.zip"

# Find Wt resources directory
WT_RES=""
if [ -n "${WT_RESOURCES_DIR}" ] && [ -d "${WT_RESOURCES_DIR}" ]; then
  WT_RES="${WT_RESOURCES_DIR}"
else
  for candidate in build_xcode build_vscode build; do
    if [ -d "${REPO_ROOT}/${candidate}/resources" ]; then
      WT_RES="${REPO_ROOT}/${candidate}/resources"
      break
    fi
  done
fi

if [ -z "${WT_RES}" ]; then
  echo "ERROR: Cannot find Wt resources directory."
  echo "  Set WT_RESOURCES_DIR env var, or build InterSpec locally first."
  exit 1
fi

echo "Using Wt resources from: ${WT_RES}"

# Stage the files that a desktop CMake build would generate into the source tree.
#
# For non-Android builds, CMake symlinks ${PROJECT_BINARY_DIR}/{data,InterSpec_resources} to the
# source dirs, so its configure_file()/deploy_js_resource() rules land in the source tree.  Android
# sets no such symlinks (SUPPORT_DIRECTORIES is empty when ANDROID), and this script runs before
# CMake anyway -- so without staging them here the assets zip silently ships without the nuclide
# database and the D3 spectrum chart, and the app is broken at runtime.
stage_generated() {
  local src="$1" dst="$2"
  if [ ! -f "${src}" ]; then
    echo "ERROR: expected generated-asset source is missing: ${src}" >&2
    echo "  (are the git submodules checked out?)" >&2
    exit 1
  fi
  if [ ! -f "${dst}" ] || [ "${src}" -nt "${dst}" ]; then
    echo "Staging ${dst#${REPO_ROOT}/}"
    mkdir -p "$(dirname "${dst}")"
    cp "${src}" "${dst}"
  fi
}

stage_generated "${REPO_ROOT}/external_libs/SandiaDecay/sandia.decay.min.xml" \
                "${REPO_ROOT}/data/sandia.decay.xml"
for _d3 in SpectrumChartD3.js SpectrumChartD3.css d3.v3.min.js; do
  stage_generated "${REPO_ROOT}/external_libs/SpecUtils/d3_resources/${_d3}" \
                  "${REPO_ROOT}/InterSpec_resources/${_d3}"
done

mkdir -p "$(dirname "${ASSETS_ZIP}")"
rm -f "${ASSETS_ZIP}"

# Zip InterSpec data from repo root
cd "${REPO_ROOT}"
zip -9 -r \
  --exclude='*OUO*' --exclude='*ouo*' \
  --exclude='*CUI*' --exclude='*cui*' --exclude='*Cui*' \
  --exclude='*.svn*' --exclude='*.git*' --exclude='*.DS_Store' \
  "${ASSETS_ZIP}" InterSpec_resources data example_spectra

# Add Wt resources.
# Normally WT_RES is a build-tree copy already slimmed by PruneWtResources.cmake.
# The font-awesome/jPlayer/bootstrap/polished excludes are belt-and-suspenders for
# the WT_RESOURCES_DIR-override path, where an un-pruned dir could be passed in
# (see target/tools/wt_resource_audit/wt_resources_manifest.json for the full, evidence-backed removal set).
cd "${WT_RES}/.."
zip -ur -9 -r \
  --exclude='resources/font-awesome/*' --exclude='resources/jPlayer/*' \
  --exclude='resources/themes/bootstrap/*' --exclude='resources/themes/polished/*' \
  --exclude='*OUO*' --exclude='*ouo*' \
  --exclude='*CUI*' --exclude='*cui*' --exclude='*Cui*' \
  --exclude='*.svn*' --exclude='*.git*' --exclude='*.DS_Store' \
  "${ASSETS_ZIP}" resources

echo "Created ${ASSETS_ZIP}"
