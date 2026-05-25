#!/bin/bash
# Builds InterSpec for Android using a Docker container.
#
# Usage:
#   ./build_android_docker.sh image           Build the SDK Docker image (once)
#   ./build_android_docker.sh build [flavor]  Build APK+AAB (default flavor: arm8Release)
#   ./build_android_docker.sh shell           Drop into the container for debugging
#
# Signing (optional):
#   export KEYSTORE_FILE=~/.android/release.keystore
#   export KEYSTORE_PASSWORD=...
#   export KEYSTORE_ALIAS_NAME=...
#   export KEYSTORE_ALIAS_PASSWORD=...
#   ./build_android_docker.sh build
#
# The Leaflet maps key can be overridden:
#   export LEAFLET_MAPS_KEY=...
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
IMAGE_NAME="interspec-android-sdk"
CACHE_DIR="${SCRIPT_DIR}/.docker-cache"
OUTPUT_DIR="${SCRIPT_DIR}/build_android_output"

WT_TAG="3.7.1"

usage() {
  echo "Usage: $0 {image|build [flavor]|shell}"
  echo ""
  echo "Commands:"
  echo "  image           Build the Docker SDK image (run once)"
  echo "  build [flavor]  Build the APK + AAB. Default flavor: arm8Release"
  echo "                  Examples: arm8Release, arm8Debug, universalRelease"
  echo "  shell           Start an interactive shell in the container"
  exit 1
}

# Build the Docker image
cmd_image() {
  echo "==> Building Docker image '${IMAGE_NAME}'..."
  docker build -t "${IMAGE_NAME}" "${SCRIPT_DIR}"
  echo "==> Docker image '${IMAGE_NAME}' built successfully."
}

# Prepare Wt resources for the assets zip
prepare_wt_resources() {
  local wt_cache="${CACHE_DIR}/wt-resources"
  if [ -d "${wt_cache}/resources" ]; then
    echo "==> Wt resources already cached."
    return
  fi

  echo "==> Fetching Wt ${WT_TAG} source for resources directory..."
  mkdir -p "${wt_cache}"
  if ! git clone --depth=1 --branch "${WT_TAG}" \
      https://github.com/emweb/wt.git "${wt_cache}/wt-src"; then
    echo "ERROR: git clone of Wt ${WT_TAG} failed -- check network/proxy configuration." >&2
    rm -rf "${wt_cache}/wt-src"
    exit 1
  fi
  if [ ! -d "${wt_cache}/wt-src/resources" ]; then
    echo "ERROR: Wt clone succeeded but resources/ is missing at ${wt_cache}/wt-src" >&2
    rm -rf "${wt_cache}/wt-src"
    exit 1
  fi
  cp -r "${wt_cache}/wt-src/resources" "${wt_cache}/resources"
  rm -rf "${wt_cache}/wt-src"
  echo "==> Wt resources cached at ${wt_cache}/resources"
}

# Run the Gradle build inside Docker
cmd_build() {
  local flavor="${1:-arm8Release}"

  # Ensure image exists
  if ! docker image inspect "${IMAGE_NAME}" > /dev/null 2>&1; then
    echo "==> Docker image not found; building it first..."
    cmd_image
  fi

  # Create persistent cache directories
  mkdir -p "${CACHE_DIR}/gradle"
  mkdir -p "${CACHE_DIR}/app-cxx"
  mkdir -p "${CACHE_DIR}/app-build"
  mkdir -p "${OUTPUT_DIR}"

  # Prepare Wt resources on host so make_assets.sh can find them
  prepare_wt_resources
  local wt_res_dir="/cache/wt-resources/resources"

  # Build the assets zip on host if it doesn't exist yet
  local assets_zip="${SCRIPT_DIR}/InterSpec/app/src/main/assets/interspec-assets.zip"
  if [ ! -f "${assets_zip}" ]; then
    echo "==> Building assets zip..."
    WT_RESOURCES_DIR="${CACHE_DIR}/wt-resources/resources" \
      "${SCRIPT_DIR}/InterSpec/make_assets.sh"
  else
    echo "==> Assets zip already exists, skipping."
  fi

  # Signing arguments
  local signing_env=()
  if [ -n "${KEYSTORE_FILE}" ] && [ -f "${KEYSTORE_FILE}" ]; then
    signing_env+=(-v "$(realpath "${KEYSTORE_FILE}"):/keystore/release.keystore:ro")
    signing_env+=(-e "KEYSTORE_LOCATION=/keystore/release.keystore")
  fi
  [ -n "${KEYSTORE_PASSWORD}" ]       && signing_env+=(-e "KEYSTORE_PASSWORD=${KEYSTORE_PASSWORD}")
  [ -n "${KEYSTORE_ALIAS_NAME}" ]     && signing_env+=(-e "KEYSTORE_ALIAS_NAME=${KEYSTORE_ALIAS_NAME}")
  [ -n "${KEYSTORE_ALIAS_PASSWORD}" ] && signing_env+=(-e "KEYSTORE_ALIAS_PASSWORD=${KEYSTORE_ALIAS_PASSWORD}")

  local leaflet_env=()
  [ -n "${LEAFLET_MAPS_KEY}" ] && leaflet_env+=(-e "LEAFLET_MAPS_KEY=${LEAFLET_MAPS_KEY}")

  # Capitalize first letter for Gradle task name (portable, works with Bash 3.x on macOS)
  local Flavor
  Flavor="$(printf '%s' "${flavor}" | awk '{print toupper(substr($0,1,1)) substr($0,2)}')"

  echo "==> Starting Gradle build: assemble${Flavor} bundle${Flavor}..."

  docker run --rm \
    -v "${REPO_ROOT}:/src:ro" \
    -v "${SCRIPT_DIR}/InterSpec:/work/InterSpec" \
    -v "${CACHE_DIR}/gradle:/root/.gradle" \
    -v "${CACHE_DIR}/app-cxx:/work/InterSpec/app/.cxx" \
    -v "${CACHE_DIR}/app-build:/work/InterSpec/app/build" \
    -v "${CACHE_DIR}/wt-resources:/cache/wt-resources:ro" \
    "${signing_env[@]}" \
    "${leaflet_env[@]}" \
    -e "WT_RESOURCES_DIR=${wt_res_dir}" \
    -w /work/InterSpec \
    "${IMAGE_NAME}" \
    bash -c "
      # Point the CMakeLists.txt to the read-only source mount
      # The gradle project's CMakeLists.txt path is relative: ../../CMakeLists.txt
      # which resolves to /work/CMakeLists.txt — we need to symlink /work's parent dirs
      ln -sfn /src/target /work/target 2>/dev/null || true
      ln -sfn /src/CMakeLists.txt /work/CMakeLists.txt 2>/dev/null || true
      ln -sfn /src/InterSpec_resources /work/InterSpec_resources 2>/dev/null || true
      ln -sfn /src/src /work/src 2>/dev/null || true
      ln -sfn /src/external_libs /work/external_libs 2>/dev/null || true
      ln -sfn /src/data /work/data 2>/dev/null || true
      ln -sfn /src/example_spectra /work/example_spectra 2>/dev/null || true

      # Create the local.properties that Gradle needs
      echo 'sdk.dir=/opt/android-sdk' > local.properties

      chmod +x gradlew
      ./gradlew assemble${Flavor} bundle${Flavor} --no-daemon
    "

  # Copy outputs
  echo "==> Copying build outputs..."
  find "${CACHE_DIR}/app-build/outputs" \( -name '*.apk' -o -name '*.aab' \) \
    -exec cp {} "${OUTPUT_DIR}/" \; 2>/dev/null || true

  echo "==> Build complete. Outputs in ${OUTPUT_DIR}/:"
  ls -lh "${OUTPUT_DIR}/"
}

# Interactive shell for debugging
cmd_shell() {
  mkdir -p "${CACHE_DIR}/gradle"
  mkdir -p "${CACHE_DIR}/app-cxx"
  mkdir -p "${CACHE_DIR}/app-build"

  docker run --rm -it \
    -v "${REPO_ROOT}:/src:ro" \
    -v "${SCRIPT_DIR}/InterSpec:/work/InterSpec" \
    -v "${CACHE_DIR}/gradle:/root/.gradle" \
    -v "${CACHE_DIR}/app-cxx:/work/InterSpec/app/.cxx" \
    -v "${CACHE_DIR}/app-build:/work/InterSpec/app/build" \
    -w /work/InterSpec \
    "${IMAGE_NAME}" \
    bash
}


case "${1}" in
  image) cmd_image ;;
  build) cmd_build "${2}" ;;
  shell) cmd_shell ;;
  *)     usage ;;
esac
