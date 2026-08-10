#!/bin/bash

# Builds the InterSpec dependency *prefix* (boost, zlib, Wt, Eigen, Ceres) inside the manylinux
# container, then verifies it and stamps it with a completion sentinel.
#
# takes three arguments: [path to InterSpec code] [dep build working dir] [prefix dir]
#  - InterSpec code is mapped to `/interspec`
#  - the dependency build working directory to `/deps_build`  (scratch; never cached)
#  - the prefix to `/deps_prefix`                             (the product; this is what gets cached)
#
# Why a prefix rather than a build directory: `actions/checkout` rewrites every source file's
# mtime, and a freshly-checked-out CMakeLists.txt forces CMake to regenerate, so a restored *build*
# tree recompiles everything anyway - measured, a warm build-directory cache compiled exactly as
# many objects as a cold one.  A prefix is a finished product with no mtime relationship to the
# source tree, so restoring it skips the work outright.  This is the same pattern the macOS and
# Windows jobs already use.
#
# The prefix MUST be built here rather than on the host: everything in it is linked into
# InterSpecAddOn.node, so it has to come from the same gcc-toolset-14 / glibc-2.28 world the app is
# compiled in.  See target/electron/linux/check_elf_compat.py.
#
# For example, to build a prefix by hand:
#  cd /tmp
#  mkdir -p build_interspec/deps_build build_interspec/deps_prefix
#  cd build_interspec
#  git clone --recursive git@github.com:sandialabs/InterSpec.git ./InterSpec_code
#  docker run --rm -it -v `pwd`/InterSpec_code:/interspec -v `pwd`/deps_build:/deps_build -v `pwd`/deps_prefix:/deps_prefix quay.io/pypa/manylinux_2_28_x86_64:latest bash -lc '/interspec/target/electron/build_linux_deps_from_docker.sh /interspec /deps_build /deps_prefix'
#  #Then compile against it with:
#  #  docker run ... build_linux_app_from_docker.sh /interspec /build_app /build_working_dir /deps_prefix
#
# Note the `bash -lc`: a login shell is what puts gcc-toolset-14 on PATH in this image.

set -euo pipefail

if [ $# -ne 3 ]; then
  echo "Error: This script requires exactly three arguments."
  echo "    [path to InterSpec code] [dep build working dir] [prefix dir]"
  exit 1
fi

InterSpecCodePath="$1"
DepsWorkDir="$2"
DepsPrefix="$3"
Sentinel="${DepsPrefix}/.interspec_deps_complete"

mkdir -p "${DepsWorkDir}" "${DepsPrefix}"


# ------------------------------------------------------------------------------------------------
# Container prerequisites.
#
# There are none to install.  quay.io/pypa/manylinux_2_28_x86_64 is AlmaLinux 8 plus the manylinux
# tooling, and dep_build_linux.sh needs only things the image already has: gcc-toolset-14 (on PATH
# for a login shell), make, cmake, git, curl, tar, gzip, awk, sha256sum, nm.
#
# So we *assert* the toolchain rather than installing it.  If cmake or the compiler were missing, a
# `yum install` would quietly hand us AlmaLinux 8's cmake 3.20 / gcc 8.5 and we would cache a
# prefix with a different ABI under the same cache key - the exact bug class the prefix migration
# exists to remove.  Failing here is far cheaper than debugging that later.
# ------------------------------------------------------------------------------------------------
_missing=""
for _cmd in cc c++ make cmake git curl tar gzip awk sha256sum nm ; do
  command -v "${_cmd}" > /dev/null 2>&1 || _missing="${_missing} ${_cmd}"
done
if [ -n "${_missing}" ]; then
  echo "Error: the container is missing:${_missing}"
  echo "  These are all expected to be present in quay.io/pypa/manylinux_2_28_x86_64."
  echo "  Do NOT paper over this with a yum install - the toolchain version is baked into the"
  echo "  cached prefix, but would not change the cache key, so the next build would silently"
  echo "  link against libraries built by a different compiler."
  exit 1
fi

echo "=== Toolchain ==="
cc --version | head -n 1
c++ --version | head -n 1
cmake --version | head -n 1

# The dependency build applies patches with `git apply` and clones Eigen/Ceres; git wants to be
#  told the (runner-owned) source tree is safe to read from as root.
git config --global --add safe.directory "${InterSpecCodePath}" || true
GIT_HASH=$(git -C "${InterSpecCodePath}" rev-parse HEAD 2>/dev/null || echo "unknown")


# ------------------------------------------------------------------------------------------------
# Build the prefix.
#
# Note dep_build_linux.sh keeps its `<name>.installed` markers in the WORKING directory, not in the
# prefix - so a restored prefix on its own does not make a re-run a no-op.  The caller has to skip
# this script entirely on a cache hit, which .github/workflows/build_app.yml does.
# ------------------------------------------------------------------------------------------------
echo "=== Building dependency prefix ==="
# Invoked via `sh` rather than executed directly, so this works when the source tree is mounted
#  read-only (which CI does) and no `chmod +x` is possible.
sh "${InterSpecCodePath}/target/patches/dep_build_linux.sh" \
  "${InterSpecCodePath}" "${DepsWorkDir}" "${DepsPrefix}"


# ------------------------------------------------------------------------------------------------
# Only now, and only because the above exited 0, decide the prefix is usable.  Everything checked
# here is something InterSpec's own configure or link step needs; catching a gap now costs seconds,
# catching it during the compile costs the whole build.
# ------------------------------------------------------------------------------------------------
echo "=== Checking the prefix ==="
_bad=""
for _f in \
    include/Wt/WObject include/Wt/Dbo/backend/Sqlite3 include/boost/version.hpp \
    include/zlib.h include/eigen3/Eigen/Core \
    lib/libwt.a lib/libwthttp.a lib/libwtdbo.a lib/libwtdbosqlite3.a \
    lib/libboost_thread.a lib/libboost_date_time.a lib/libboost_system.a \
    lib/libboost_filesystem.a lib/libboost_program_options.a lib/libboost_regex.a \
    lib/libboost_random.a lib/libz.a lib/libceres.a \
    share/Wt/resources/form.css ; do
  [ -e "${DepsPrefix}/${_f}" ] || _bad="${_bad} ${_f}"
done
if [ -n "${_bad}" ]; then
  echo "Error: the prefix is missing:${_bad}"
  exit 1
fi

# The CMake package files find_package() actually consumes.  b2 writes Boost's under a
#  version-stamped directory, so glob rather than naming a path.
ls -d "${DepsPrefix}"/lib/cmake/Boost-* > /dev/null 2>&1 \
  || { echo "Error: no BoostConfig.cmake under ${DepsPrefix}/lib/cmake/Boost-*"; exit 1; }

# `lib`, not `lib64`: dep_build_linux.sh passes -DCMAKE_INSTALL_LIBDIR=lib to every CMake
#  sub-build, because GNUInstallDirs would otherwise pick lib64 on this (RHEL-family) image while
#  boost and Wt install into lib.  A split prefix mostly works but breaks the two places that
#  hardcode `lib` - the libz.so cleanup below, and cmake/FindWt.cmake's `HINTS <inc>/../lib`.
#  Check here so a layout change is a clear error rather than a mysterious find_package failure.
for _cfg in lib/cmake/Ceres/CeresConfig.cmake share/eigen3/cmake/Eigen3Config.cmake ; do
  if [ ! -e "${DepsPrefix}/${_cfg}" ]; then
    echo "Error: ${DepsPrefix}/${_cfg} is missing."
    echo "  If it turned up under lib64/ instead, -DCMAKE_INSTALL_LIBDIR=lib did not take effect."
    find "${DepsPrefix}" -name 'CeresConfig.cmake' -o -name 'Eigen3Config.cmake' 2>/dev/null || true
    exit 1
  fi
done

# A shared libz in the prefix would show up as `DT_NEEDED libz.so.1` on InterSpecAddOn.node, which
#  check_elf_compat.py rejects on purpose (it is not on the allowlist).  dep_build_linux.sh deletes
#  them after installing zlib; make sure it actually did.
if ls "${DepsPrefix}"/lib/libz.so* > /dev/null 2>&1; then
  echo "Error: ${DepsPrefix}/lib contains a shared libz; it must be static-only."
  ls -la "${DepsPrefix}"/lib/libz.so*
  exit 1
fi

# cmake/FindWt.cmake hands InterSpec bare .a paths with no link interface, so anything Wt itself
#  linked against has to already be *inside* the archive.  Wt falls back to its bundled sqlite3
#  amalgamation when no system sqlite3 is found, which is the case in this image.  If that ever
#  changes, InterSpec's link fails with undefined sqlite3_* references, and this is the line that
#  explains why.
#  NB: `grep -c`, not `grep -q`.  `grep -q` exits the moment it matches, which closes the pipe and
#  kills `nm` (megabytes of output) with SIGPIPE - and under `set -o pipefail` that 141 becomes the
#  pipeline's status, so a successful match reads as a failure.  `grep -c` consumes all of nm's
#  output, and the `|| true` is for the genuine no-match case, where grep exits 1.
_sqlite_syms=$(nm -g --defined-only "${DepsPrefix}/lib/libwtdbosqlite3.a" 2>/dev/null | grep -c " T sqlite3_open" || true)
if [ "${_sqlite_syms:-0}" -eq 0 ]; then
  echo "WARNING: libwtdbosqlite3.a does not define sqlite3_open - Wt appears to have linked against"
  echo "         a system sqlite3.  cmake/FindWt.cmake will not propagate that, so expect undefined"
  echo "         sqlite3_* symbols when InterSpecAddOn.node is linked."
else
  echo "Wt's sqlite3 is bundled into libwtdbosqlite3.a (${_sqlite_syms} sqlite3_open* symbols) - good."
fi

# The sentinel is both a "this prefix is complete" flag and a provenance record.  build_app.yml
#  prints it on every run - including cache hits - so "which prefix am I linking against?" is
#  answerable from the run log rather than by inference.
{
  echo "interspec_git_hash=${GIT_HASH}"
  echo "built_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "prefix_path=${DepsPrefix}"
  echo "arch=$(uname -m)"
  echo "cc=$(cc --version | head -n 1)"
  echo "cmake=$(cmake --version | head -n 1)"
  echo "dep_script_sha256=$(sha256sum "${InterSpecCodePath}/target/patches/dep_build_linux.sh" | awk '{print $1}')"
  echo "wt_patch_sha256=$(sha256sum "${InterSpecCodePath}/target/patches/wt/3.7.1/NormalBuild/wt_3.7.1_git.patch" | awk '{print $1}')"
} > "${Sentinel}"

echo "=== Prefix complete ==="
cat "${Sentinel}"
du -sh "${DepsPrefix}"
