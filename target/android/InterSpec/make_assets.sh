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

rm -f "${ASSETS_ZIP}"

# Zip InterSpec data from repo root
cd "${REPO_ROOT}"
zip -9 -r \
  --exclude='*OUO*' --exclude='*ouo*' \
  --exclude='*CUI*' --exclude='*cui*' --exclude='*Cui*' \
  --exclude='*.svn*' --exclude='*.git*' --exclude='*.DS_Store' \
  "${ASSETS_ZIP}" InterSpec_resources data example_spectra

# Add Wt resources
cd "${WT_RES}/.."
zip -ur -9 -r \
  --exclude='*OUO*' --exclude='*ouo*' \
  --exclude='*CUI*' --exclude='*cui*' --exclude='*Cui*' \
  --exclude='*.svn*' --exclude='*.git*' --exclude='*.DS_Store' \
  "${ASSETS_ZIP}" resources

echo "Created ${ASSETS_ZIP}"
