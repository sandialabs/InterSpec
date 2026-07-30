#!/usr/bin/env python3
"""Assert the Linux binaries we ship stay within the runtime floor we advertise.

This exists because of https://github.com/sandialabs/InterSpec/issues/51.  Every Electron 41
Linux build segfaulted in `dl_init` before `main()`, because `InterSpecAddOn.node` linked
libstdc++ statically: libstdc++'s own locale-facet initialization became just another static
initializer inside that shared object, and the toolchain ordered it after Wt's, so the
namespace-scope `boost::basic_regex` in Wt's CgiParser.C called through `std::collate<char>`
and landed in the wrong vtable slot.  Linking libstdc++ dynamically makes it a real shared
object whose `dl_init` the loader is required to run first.

So the two things worth checking on every build are (a) that libstdc++ really is dynamic, and
(b) that going dynamic did not quietly raise the floor above what the release notes promise.
The release build uses quay.io/pypa/manylinux_2_28_x86_64, whose compiler is Red Hat
gcc-toolset-14; Developer Toolset links the *base* AlmaLinux 8 runtime and pulls only newer
ABI bits from libstdc++_nonshared.a, so the ceilings below are what that base provides.

Usage:  check_elf_compat.py <directory-or-file> [...]

Only binaries InterSpec builds are checked.  Electron's own prebuilt binaries are reported for
information but never failed on - we do not control how they were linked.
"""

import argparse
import os
import re
import subprocess
import sys


# What the AlmaLinux 8 base runtime provides; equivalently, the floor advertised on the release
# page (glibc >= 2.28: RHEL 8, Ubuntu 20.04+, Debian 11+, Fedora 30+).
MAX_VERSIONS = {
  "GLIBC": (2, 28),
  "GLIBCXX": (3, 4, 25),
  "CXXABI": (1, 3, 11),
  "GCC": (7, 0, 0),
}

# Shared libraries InterSpec's own binaries may depend on.  Anything else means a dependency
# crept in that end-user systems are not guaranteed to have.
ALLOWED_NEEDED = {
  "libstdc++.so.6",
  "libgcc_s.so.1",
  "libc.so.6",
  "libm.so.6",
  "libdl.so.2",
  "librt.so.1",
  "libpthread.so.0",
  # The dynamic loader itself; shows up as NEEDED for some link configurations and can never be
  # a portability problem, since a dynamically-linked binary cannot run without it.
  "ld-linux-x86-64.so.2",
}

# Note libz.so.1 is deliberately NOT in the list above.  It used to appear in DT_NEEDED as pure
# over-link: zlib is compiled from source by FetchContent and linked statically (the addon
# defines deflate/inflate/zlibVersion itself and imports no zlib symbol), but Wt ran its own
# find_package(ZLIB) against /usr - before our static target existed - and put the system shared
# library into wthttp's link interface.  Setting HTTP_WITH_ZLIB=OFF in
# cmake/FetchInterSpecDeps.cmake stops that at the source.  If libz.so.1 reappears, something is
# genuinely linking the system zlib and that deserves a look rather than an allowlist entry.

# Must be dynamically linked, not absorbed into the binary - this is the issue #51 guard.
REQUIRED_NEEDED = { "libstdc++.so.6" }

# Binaries InterSpec builds, as opposed to the Electron runtime we merely redistribute.
OURS = re.compile( r"(^InterSpecAddOn\.node$)|(^libInterSpec\.so)|(^InterSpec_batch$)" )

VERSION_RE = re.compile( r"\b(GLIBC|GLIBCXX|CXXABI|GCC)_(\d+(?:\.\d+)*)\b" )
NEEDED_RE = re.compile( r"\(NEEDED\).*\[([^\]]+)\]" )


def readelf( args, path ):
  """Return readelf output, or None if the file is not an ELF object."""
  try:
    res = subprocess.run( [ "readelf" ] + args + [ path ],
                          capture_output=True, text=True, check=False )
  except FileNotFoundError:
    sys.exit( "error: `readelf` not found - install binutils" )

  if res.returncode != 0 or "Not an ELF file" in res.stderr:
    return None
  return res.stdout


def is_elf( path ):
  try:
    with open( path, "rb" ) as f:
      return f.read( 4 ) == b"\x7fELF"
  except OSError:
    return False


def required_versions( path ):
  """Max required version per prefix, from the .gnu.version_r (needed) entries only.

  `readelf -V` prints both what the object *provides* (.gnu.version_d) and what it *requires*
  (.gnu.version_r).  Only the latter constrains where the binary can run, so parse just the
  "Version needs section" blocks.
  """
  out = readelf( [ "-V", "-W" ], path )
  if out is None:
    return {}

  maxima = {}
  in_needs = False
  for line in out.splitlines():
    if line.startswith( "Version needs section" ):
      in_needs = True
      continue
    if line.startswith( "Version definition section" ) or line.startswith( "Version symbols section" ):
      in_needs = False
      continue
    if not in_needs:
      continue

    for prefix, ver in VERSION_RE.findall( line ):
      parsed = tuple( int(p) for p in ver.split(".") )
      if parsed > maxima.get( prefix, (0,) ):
        maxima[prefix] = parsed

  return maxima


def needed_libs( path ):
  out = readelf( [ "-d", "-W" ], path )
  if out is None:
    return set()
  return set( NEEDED_RE.findall( out ) )


def fmt( prefix, ver ):
  return "%s_%s" % (prefix, ".".join( str(p) for p in ver ))


def check( path, strict ):
  """Return a list of problem strings; empty means the binary is fine."""
  problems = []
  name = os.path.basename( path )

  maxima = required_versions( path )
  needed = needed_libs( path )

  detail = ", ".join( sorted( fmt(p, v) for p, v in maxima.items() ) ) or "none"
  print( "  %-28s requires %s" % (name, detail) )
  print( "  %-28s NEEDED   %s" % ("", ", ".join( sorted( needed ) ) or "none") )

  if not strict:
    return problems

  for prefix, ver in sorted( maxima.items() ):
    ceiling = MAX_VERSIONS.get( prefix )
    if ceiling and ver > ceiling:
      problems.append( "%s requires %s but the floor allows at most %s"
                       % (name, fmt(prefix, ver), fmt(prefix, ceiling)) )

  for lib in sorted( REQUIRED_NEEDED - needed ):
    problems.append( "%s does not link %s dynamically - this is the issue #51 regression;"
                     " see cmake/CxxRuntimePolicy.cmake" % (name, lib) )

  for lib in sorted( needed - ALLOWED_NEEDED ):
    problems.append( "%s depends on %s, which is not in the allowlist" % (name, lib) )

  return problems


def main():
  parser = argparse.ArgumentParser( description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter )
  parser.add_argument( "paths", nargs="+", help="Files, or directories to search recursively" )
  args = parser.parse_args()

  candidates = []
  for path in args.paths:
    if os.path.isdir( path ):
      for root, _dirs, files in os.walk( path ):
        for f in sorted( files ):
          full = os.path.join( root, f )
          if not os.path.islink( full ) and is_elf( full ):
            candidates.append( full )
    elif is_elf( path ):
      candidates.append( path )
    else:
      sys.exit( "error: %s is not an ELF file or directory" % path )

  if not candidates:
    sys.exit( "error: no ELF binaries found in %s" % ", ".join( args.paths ) )

  ours = [ p for p in candidates if OURS.match( os.path.basename( p ) ) ]
  theirs = [ p for p in candidates if p not in ours ]

  if not ours:
    sys.exit( "error: none of the %d ELF binaries found look like InterSpec's own"
              " (expected InterSpecAddOn.node); wrong directory?" % len(candidates) )

  problems = []

  print( "InterSpec binaries (enforced):" )
  for path in ours:
    problems += check( path, strict=True )

  if theirs:
    print( "\nRedistributed Electron binaries (informational only):" )
    for path in theirs:
      check( path, strict=False )

  print()
  if problems:
    print( "FAIL: %d problem(s)" % len(problems) )
    for p in problems:
      print( "  - " + p )
    return 1

  print( "OK: %d InterSpec binary(s) within the declared runtime floor" % len(ours) )
  return 0


if __name__ == "__main__":
  sys.exit( main() )
