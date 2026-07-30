# How InterSpec links the C++ runtime (libstdc++ / libgcc) on Linux.
#
# This is deliberately the single place the policy is expressed.  It used to live as a
# `--CDCMAKE_SHARED_LINKER_FLAGS="-static-libgcc -static-libstdc++"` string inside
# target/electron/build_linux_app_from_docker.sh, which was invisible to CMake, applied only
# to SHARED link lines (not EXE or MODULE), and left every other Linux consumer - the
# target/docker container builds, the plain BUILD_AS_LOCAL_SERVER build, a distro packager -
# on a different policy.
#
# The default is now to link libstdc++ DYNAMICALLY, which fixes
# https://github.com/sandialabs/InterSpec/issues/51: every Electron 41 Linux build segfaulted
# before main().  With libstdc++ linked statically into InterSpecAddOn.node, libstdc++'s own
# locale-facet initialization is just another static initializer in that same shared object,
# and the toolchain ordered it *after* Wt's.  So the namespace-scope boost::basic_regex in
# Wt's CgiParser.C ran first, called through std::collate<char>, and landed in the wrong
# vtable slot (std::codecvt<char16_t,...>::do_unshift).  Linking libstdc++ dynamically makes
# it a real shared object, so the loader guarantees its dl_init runs before ours.
#
# Dynamic linking does not raise the supported-OS floor.  The release build uses
# quay.io/pypa/manylinux_2_28_x86_64, whose compiler is Red Hat gcc-toolset-14; Developer
# Toolset links the *base* AlmaLinux 8 libstdc++.so.6 and pulls only newer ABI bits from
# libstdc++_nonshared.a, so the produced binaries require at most GLIBCXX_3.4.25 /
# CXXABI_1.3.11 / GLIBC_2.28 - satisfied by every distro we claim (RHEL 8, Ubuntu 20.04+,
# Debian 11+, Fedora 30+).  target/electron/linux/check_elf_compat.py asserts this in CI, so
# a regression here fails the build instead of shipping a segfault.
#
# If a future toolchain ever forces static linking back on, the initialization-order hazard
# returns and Wt needs patching as well: make the CgiParser.C regex a function-local static
# (add it to target/patches/wt/3.7.1/FetchContent/wt_3.7.1_FetchContent.git.patch).

option( INTERSPEC_STATIC_CXX_RUNTIME
        "Statically link libstdc++/libgcc into InterSpec (see comment in cmake/CxxRuntimePolicy.cmake)"
        OFF )

# A fully-static build supplies its own runtime wholesale via CMAKE_EXE_LINKER_FLAGS, so it
# does not need - and must not duplicate - these flags.
if( CONTAINER_MUSL )
  set( INTERSPEC_STATIC_CXX_RUNTIME OFF )
endif( CONTAINER_MUSL )

if( INTERSPEC_STATIC_CXX_RUNTIME
    AND UNIX AND NOT APPLE AND NOT ANDROID
    AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" )
  # PUBLIC on the static InterSpecLib is what makes this consistent: the options propagate
  # through INTERFACE_LINK_OPTIONS to every consumer - InterSpecAddOn.node, InterSpecExe,
  # InterSpec_batch, the unit tests - whatever their link type.  add_link_options() would not
  # work here, because target/electron/CMakeLists.txt is the *parent* scope of this project
  # (it does add_subdirectory(../..)), so directory-scoped options set here would never reach
  # the InterSpecAddOn target.
  target_link_options( InterSpecLib PUBLIC -static-libgcc -static-libstdc++ )
  message( STATUS "InterSpec: statically linking libstdc++/libgcc" )
endif()
