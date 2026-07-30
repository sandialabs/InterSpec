#!/bin/sh
# Launcher for the extracted-tarball build of InterSpec.
#
# The .deb and .rpm do not use this: their postinst makes `chrome-sandbox` setuid-root and
# installs an AppArmor profile, so the real binary can be run directly.  A tarball has no
# install step, so nothing can grant those privileges, and Chromium aborts outright rather than
# quietly running unsandboxed.  Rather than hard-code `--no-sandbox` and give up the sandbox even
# where it would work, probe for the two things it needs:
#
#   1. chrome-sandbox owned by root with the setuid bit - tar preserves the bit, but only if the
#      archive was unpacked by root, which is not the normal case.
#   2. Unprivileged user namespaces - blocked by default on Ubuntu 23.10+ for binaries with no
#      AppArmor profile (kernel.apparmor_restrict_unprivileged_userns=1).
#
# See https://github.com/sandialabs/InterSpec/issues/51 for the failure modes this avoids.

set -e

# Resolve the real directory of this script, following symlinks, so the tarball can be extracted
# anywhere and still be launched through a symlink on $PATH.
SELF="$0"
while [ -L "$SELF" ]; do
  LINK=$( ls -ld -- "$SELF" | sed 's/.* -> //' )
  case "$LINK" in
    /*) SELF="$LINK" ;;
    *)  SELF="$( dirname -- "$SELF" )/$LINK" ;;
  esac
done
DIR=$( cd -- "$( dirname -- "$SELF" )" && pwd -P )

BIN="$DIR/interspec-bin"

if [ ! -x "$BIN" ]; then
  echo "InterSpec: cannot find the executable at $BIN" >&2
  exit 1
fi

# If the user asked for a specific sandbox setting, respect it and get out of the way.
for arg in "$@"; do
  case "$arg" in
    --no-sandbox|--no-sandbox=*|--enable-sandbox)
      exec "$BIN" "$@"
      ;;
  esac
done

# The two sandbox mechanisms are alternatives, not both required - the same logic
# electron-builder's own .deb/.rpm post-install script uses to decide whether chrome-sandbox
# needs the setuid bit at all.
sandbox_usable() {
  # Preferred: the namespace sandbox, which needs no setuid helper.  Blocked on Ubuntu 23.10+
  # for binaries with no AppArmor profile, which is exactly the tarball's situation.
  if [ -e /proc/self/ns/user ] \
     && command -v unshare >/dev/null 2>&1 \
     && unshare -Ur true >/dev/null 2>&1; then
    return 0
  fi

  # Fallback: the SUID sandbox, which Chromium uses when user namespaces are unavailable.
  if [ -u "$DIR/chrome-sandbox" ] \
     && [ "$( stat -c %u "$DIR/chrome-sandbox" 2>/dev/null )" = "0" ]; then
    return 0
  fi

  return 1
}

if sandbox_usable; then
  exec "$BIN" "$@"
fi

echo "InterSpec: no Chromium sandbox is available here - unprivileged user namespaces are" >&2
echo "  restricted and chrome-sandbox is not setuid-root - so starting with --no-sandbox." >&2
echo "  Installing the .deb or .rpm package instead lets InterSpec run sandboxed." >&2

exec "$BIN" --no-sandbox "$@"
