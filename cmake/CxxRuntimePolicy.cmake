# How InterSpec links the C++ runtime (libstdc++ / libgcc) on Linux.
#
# This is deliberately the single place the policy is expressed.  It used to live as a
# `--CDCMAKE_SHARED_LINKER_FLAGS="-static-libgcc -static-libstdc++"` string inside
# target/electron/build_linux_app_from_docker.sh, which applied only to SHARED link lines (not
# EXE or MODULE) and left every other Linux consumer - the target/docker container builds, the
# plain BUILD_AS_LOCAL_SERVER build, a distro packager - on a different policy.
#
# Note `--CD<VAR>` becomes a `-D<VAR>` *cache* entry, so simply dropping the flag from the build
# script does not undo it in an existing build directory.  CI restores its build directory from
# a cache, so the first run after that change reconfigured on top of a CMakeCache.txt that still
# had the flag, statically linked libstdc++ again, and reproduced issue #51 - caught only by
# check_elf_compat.py.  Hence the sanitizing block below: this file has to *enforce* the policy,
# not merely express a default.
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

if( UNIX AND NOT APPLE AND NOT ANDROID
    AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang" )

  if( CONTAINER_MUSL )
    # A fully-static build owns its runtime linkage wholesale via CMAKE_EXE_LINKER_FLAGS (set in
    # the top-level CMakeLists.txt), so leave it entirely alone - both the option below and the
    # sanitizing, which would otherwise strip flags that build genuinely wants.
    message( STATUS "InterSpec: CONTAINER_MUSL - leaving C++ runtime linkage to the static build" )
  elseif( INTERSPEC_STATIC_CXX_RUNTIME )
    # PUBLIC on the static InterSpecLib is what makes this consistent: the options propagate
    # through INTERFACE_LINK_OPTIONS to every consumer - InterSpecAddOn.node, InterSpecExe,
    # InterSpec_batch, the unit tests - whatever their link type.  add_link_options() would not
    # work here, because target/electron/CMakeLists.txt is the *parent* scope of this project
    # (it does add_subdirectory(../..)), so directory-scoped options set here would never reach
    # the InterSpecAddOn target.
    target_link_options( InterSpecLib PUBLIC -static-libgcc -static-libstdc++ )
    message( STATUS "InterSpec: statically linking libstdc++/libgcc" )
  else()
    # Actively strip the flags rather than just not adding them, so an inherited value - a
    # restored CI cache, a stale local build dir, a packager's environment - cannot silently
    # re-introduce issue #51.  Warn loudly, because it means something upstream disagrees with
    # this policy and we are overriding it.
    foreach( _flag_var IN ITEMS CMAKE_EXE_LINKER_FLAGS CMAKE_SHARED_LINKER_FLAGS
                                CMAKE_MODULE_LINKER_FLAGS CMAKE_CXX_FLAGS )
      if( "${${_flag_var}}" MATCHES "-static-lib(stdc\\+\\+|gcc)" )
        string( REGEX REPLACE "-static-libstdc\\+\\+" "" _cleaned "${${_flag_var}}" )
        string( REGEX REPLACE "-static-libgcc" "" _cleaned "${_cleaned}" )
        string( STRIP "${_cleaned}" _cleaned )

        get_property( _doc CACHE ${_flag_var} PROPERTY HELPSTRING )
        set( ${_flag_var} "${_cleaned}" CACHE STRING "${_doc}" FORCE )
        set( ${_flag_var} "${_cleaned}" )

        message( WARNING
          "InterSpec: removed -static-libstdc++/-static-libgcc from ${_flag_var}."
          "  Statically linking libstdc++ into the Electron addon causes the startup crash in"
          "  issue #51.  Set INTERSPEC_STATIC_CXX_RUNTIME=ON if you really want it." )
      endif()
    endforeach()
  endif()

endif()
