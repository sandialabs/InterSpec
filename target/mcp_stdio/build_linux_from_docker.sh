#!/bin/bash

# Builds the stdio MCP server (target/mcp_stdio) for Linux, inside the manylinux container, against a
#  dependency prefix built by target/electron/build_linux_deps_from_docker.sh, and installs it.
#
# takes four arguments: [path to InterSpec code] [CMake build dir] [install dir] [dependency prefix]
#  - InterSpec code is mapped to `/interspec` (read-write: the InterSpec build writes into its data/ and
#    InterSpec_resources/ directories, through symlinks from the build directory)
#  - the build directory to `/build_stdio`
#  - the install directory to `/stdio_install`; the server is installed to <install dir>/InterSpec-mcp-stdio
#  - the dependency prefix to `/deps_prefix` (read-only is fine)
#
# Building in quay.io/pypa/manylinux_2_28_x86_64 is what makes the executable run on any distro with
#  glibc 2.28 or newer (RHEL 8, Debian 10, Ubuntu 20.04, and later) - the same floor as the Electron
#  app.  It links libstdc++ and libgcc statically, so glibc is all it needs from the system; libcurl
#  (for the assistant's LLM calls only) is loaded at run time if it is there.
#
# For example, to build by hand (after building a prefix as described in
#  build_linux_deps_from_docker.sh), where <image> must be the image the prefix was built in - the
#  prefix is ABI-bound to its gcc-toolset; CI uses the default of .github/actions/linux_deps_prefix:
#  docker run --rm -v `pwd`/InterSpec_code:/interspec -v `pwd`/build_stdio:/build_stdio -v `pwd`/stdio_install:/stdio_install -v `pwd`/deps_prefix:/deps_prefix:ro <image> bash -lc '/interspec/target/mcp_stdio/build_linux_from_docker.sh /interspec /build_stdio /stdio_install /deps_prefix'
#
# Note the `bash -lc`: a login shell is what puts gcc-toolset-14 on PATH in this image.

set -euo pipefail

if [ $# -ne 4 ]; then
  echo "Error: This script requires exactly four arguments."
  echo "    [path to InterSpec code] [CMake build dir] [install dir] [dependency prefix]"
  exit 1
fi

InterSpecCodePath="$1"
BuildDir="$2"
InstallDir="$3"
DepsPrefix="$4"

if [ ! -f "${DepsPrefix}/.interspec_deps_complete" ]; then
  echo "Error: '${DepsPrefix}' is not a completed InterSpec dependency prefix (no"
  echo "  .interspec_deps_complete sentinel); build one with"
  echo "  target/electron/build_linux_deps_from_docker.sh."
  exit 1
fi
echo "Using dependency prefix ${DepsPrefix}:"
cat "${DepsPrefix}/.interspec_deps_complete"

echo "=== Toolchain ==="
c++ --version | head -n 1
cmake --version | head -n 1

git config --global --add safe.directory "${InterSpecCodePath}" || true

# Headers only: the LLM transport dlopen()s libcurl.so.4 at run time, so nothing links it.  Without
#  them CMake would only warn, and fall back to Wt::Http::Client - which has no HTTPS, as the prefix's
#  Wt is built without SSL.  (No `yum update`: it could upgrade gcc-toolset past the compiler the
#  prefix was built with.)
yum install -y libcurl-devel

. "${InterSpecCodePath}/target/dep_build/manylinux_install_ccache.sh"
manylinux_install_ccache

LAUNCHER_ARGS=()
if command -v ccache > /dev/null 2>&1; then
  ccache --version | head -n 1
  ccache --zero-stats > /dev/null 2>&1 || true
  LAUNCHER_ARGS=( -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache )
fi

# Unix Makefiles, as the image has no ninja.  ZLIB_ROOT keeps CMake from picking up the image's
#  /usr/lib64/libz.so over the prefix's static one.
cmake -S "${InterSpecCodePath}/target/mcp_stdio" -B "${BuildDir}" -G "Unix Makefiles" \
  -DCMAKE_BUILD_TYPE=Release \
  -DInterSpec_FETCH_DEPENDENCIES=OFF \
  -DCMAKE_PREFIX_PATH="${DepsPrefix}" \
  -DZLIB_ROOT="${DepsPrefix}" \
  -DCURL_INCLUDE_DIR=/usr/include \
  ${LAUNCHER_ARGS[@]+"${LAUNCHER_ARGS[@]}"}

cmake --build "${BuildDir}" --target InterSpecMcpStdio test_quickjs_sandbox --parallel "$(nproc)"

# The sandbox the assistant's run_javascript runs in; a security boundary, so worth the one link.
ctest --test-dir "${BuildDir}" -R quickjs_sandbox --output-on-failure

rm -rf "${InstallDir}/InterSpec-mcp-stdio"
cmake --install "${BuildDir}" --component McpStdio --prefix "${InstallDir}/InterSpec-mcp-stdio"

Exe="${InstallDir}/InterSpec-mcp-stdio/bin/InterSpec-mcp-stdio"

# glibc 2.28 ceiling, and a static C++ runtime.
python3 "${InterSpecCodePath}/target/electron/linux/check_elf_compat.py" --static-cxx-runtime "${Exe}"

# The libcurl backend names the libraries it loads; without it, LLM calls would have no HTTPS.  (Look
#  for the fallback name: "libcurl.so.4" alone also appears in an error message of mcp_stdio_main.cpp.)
if ! strings "${Exe}" | grep -c "libcurl-gnutls.so.4" > /dev/null; then
  echo "Error: ${Exe} has no libcurl HTTP backend (src/NativeHttpClient_curl.cpp was not built)."
  exit 1
fi

if command -v ccache > /dev/null 2>&1; then
  echo "=== ccache statistics for this build ==="
  ccache --show-stats --verbose 2>/dev/null || ccache --show-stats || true
fi

echo "Done - the server is in ${InstallDir}/InterSpec-mcp-stdio"
