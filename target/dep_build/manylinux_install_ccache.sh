# Sourced (not executed) by the scripts that compile InterSpec inside the manylinux container:
#  target/electron/build_linux_app_from_docker.sh and target/mcp_stdio/build_linux_from_docker.sh.
#
# Defines manylinux_install_ccache, which puts ccache in /usr/local/bin if it is not already on PATH.
#  Unlike a cached build directory, ccache is content-addressed, so it is immune to `actions/checkout`
#  rewriting every source mtime - which is the reason caching the build tree never saved CI any
#  compilation at all.
#
# ccache is not in the AlmaLinux 8 repos this image uses (el8 moved it to EPEL, which manylinux does
#  not enable), so take the official release binary: the musl-static build, which is fully static and
#  so has no glibc floor of its own, as a .tar.gz so we need no `xz` either.  The version and checksums
#  are pinned; bump them deliberately.
#
# Never fails: ccache is an optimization, not a correctness requirement.  Callers check
#  `command -v ccache` afterwards to decide whether to use it.  Everything is written as `if ... && ...`
#  chains, so it is safe under both `set -e` and an ERR trap.
manylinux_install_ccache()
{
  local _cc_version="4.13.6"
  local _cc_arch=""
  local _cc_sha256=""

  if command -v ccache > /dev/null 2>&1; then
    return 0
  fi

  case "$(uname -m)" in
    x86_64)  _cc_arch="x86_64";  _cc_sha256="09e0547a0c3b250a76675c33130366f1399f3580842fb360c052520d56214ead" ;;
    aarch64) _cc_arch="aarch64"; _cc_sha256="bff0e0c19165db8627c85c36b0885b3b180659eda67a028298d26565aad52f56" ;;
    *)       _cc_arch="" ;;
  esac

  if [ -z "${_cc_arch}" ]; then
    echo "WARNING: no pinned ccache build for $(uname -m) - compiling without ccache"
    return 0
  fi

  local _cc_dir="ccache-${_cc_version}-linux-${_cc_arch}-musl-static"
  local _cc_tar="/tmp/${_cc_dir}.tar.gz"
  local _cc_url="https://github.com/ccache/ccache/releases/download/v${_cc_version}/${_cc_dir}.tar.gz"
  echo "Installing ccache ${_cc_version} (${_cc_arch})"
  if curl -fsSL "${_cc_url}" -o "${_cc_tar}" \
     && echo "${_cc_sha256}  ${_cc_tar}" | sha256sum --check --status \
     && tar -xzf "${_cc_tar}" -C /tmp \
     && install -m 0755 "/tmp/${_cc_dir}/ccache" /usr/local/bin/ccache ; then
    echo "Installed ccache from the official release tarball"
  elif yum install -y epel-release > /dev/null 2>&1 && yum install -y ccache > /dev/null 2>&1 ; then
    echo "Installed ccache from EPEL"
  else
    echo "WARNING: could not install ccache - compiling without it"
  fi
  return 0
}
