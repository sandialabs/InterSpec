cmake_minimum_required(VERSION 3.20...3.31)

set( SpecUtils_JAVA_SWIG OFF CACHE INTERNAL "" )
set( SpecUtils_D3_SUPPORT_FILE_STATIC OFF CACHE INTERNAL "" )

set( ANDROID_CPP_FEATURES "rtti exceptions" CACHE INTERNAL "" )
set( ANDROID_PIE ON CACHE INTERNAL "" )
set( ANDROID_STL "c++_static" CACHE INTERNAL "" )
set( BOOST_EXCLUDE_LIBRARIES "locale" CACHE INTERNAL "" )


set( BOOST_LOCALE_ENABLE_ICONV OFF CACHE INTERNAL "" )
set( BOOST_LOCALE_ENABLE_ICU OFF CACHE INTERNAL "" )

set( BUILD_TESTING OFF CACHE INTERNAL "" )

set( InterSpec_FETCH_DEPENDENCIES ON CACHE INTERNAL "" )


set( BUILD_AS_LOCAL_SERVER OFF CACHE INTERNAL "" )
set( TRY_TO_STATIC_LINK ON CACHE INTERNAL "" )
set( MAX_SPECTRUM_MEMMORY_SIZE_MB 32 CACHE INTERNAL "" )
set( USE_DB_TO_STORE_SPECTRA ON CACHE INTERNAL "" )
set( USE_SPECRUM_FILE_QUERY_WIDGET OFF CACHE INTERNAL "" )
set( USE_REL_ACT_TOOL ON CACHE INTERNAL "" )
set( USE_LEAFLET_MAP ON CACHE INTERNAL "" )


message( "ANDROID_NDK_SYSROOT=${ANDROID_NDK_SYSROOT}")

# Wt 4.13.2's src/http/Android.C guards its body (including the
# preventRemoveOfSymbolsDuringLinking() definition that WServer.C
# references under ANDROID) with `#ifdef WT_THREADED`, but never
# includes WConfig.h before that check.  Defining WT_THREADED at the
# command line makes the guard visible and lets the file compile.
add_compile_definitions( WT_THREADED )

project( InterSpecAndroid VERSION 1.0.14 )

# NDK r28 with minSdk 21 enables NEON by default for arm64-v8a.
# The __ANDROID__ macro is defined by the NDK toolchain, so no need to add it manually.

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/LibInterSpec )

# TODO: test if ${log-lib} can replace log below
add_library( InterSpecAppLib SHARED android.cpp )
target_link_libraries( InterSpecAppLib PUBLIC Boost::filesystem wt InterSpecLib log )
target_include_directories( InterSpecAppLib PUBLIC ${wt_BINARY_DIR} )

# Wt 4 bundles its own copy of Nayuki's qrcodegen (for Auth::Mfa::TotpQrCode)
# inside libwt.a, and InterSpec links its own copy from
# external_libs/QR-Code-generator/.  Both are byte-compatible at the
# implementation level (same upstream source), but lld's strict default
# refuses the duplicate symbols when linking the InterSpecAppLib .so.
# Desktop builds (Mach-O / dynamic linking) tolerate this automatically;
# on Android ELF + static linking we have to tell lld to accept the
# multiple definitions and pick the first.
target_link_options( InterSpecAppLib PRIVATE -Wl,--allow-multiple-definition )
