find_package(Git)
include(FetchContent)

# ------------------------------------------------------------------
# Prefer local copies that CMake has already unpacked to build/_deps
# ------------------------------------------------------------------
macro(use_local_dep dep_name)
  set(_local_dir "${CMAKE_BINARY_DIR}/_deps/${dep_name}-src")
  if(EXISTS "${_local_dir}/CMakeLists.txt")
    message(STATUS "Using cached ${dep_name} from ${_local_dir}")
    set(FETCHCONTENT_SOURCE_DIR_${dep_name} "${_local_dir}" CACHE PATH "" FORCE)
    set(FETCHCONTENT_UPDATES_DISCONNECTED_${dep_name} ON CACHE BOOL "" FORCE)
  endif()
endmacro()

foreach(_dep IN ITEMS boost wt zlib eigen ceres-solver wxWidgets)
  use_local_dep(${_dep})
endforeach()

# We'll set the install directory to somewhere not /usr/local or something, so we dont
#  mess the system filesystem up; Wt and zlib will install 
#  If you want to install the packages, it looks like doing the following command
#   should place them in prefix directory - although I think I patched the install 
#   functionality out of zlib and Wt... I should restore.
#     cmake --build . --config Release -j8 --target install
set( CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/prefix" CACHE INTERNAL "" )

set(BOOST_ENABLE_CMAKE ON CACHE INTERNAL "" )


# Currently installing dependencies to "prefix" of the build directory is a work in progress, 
#  so just doing it for iOS build, at the moment.... Boost and Wt seem to be not working - not sure if this is something I did, or not
if( InterSpec_IOS )
  set( INSTALL_DEPENDENCIES_IN_BUILD_DIR ON )
else( InterSpec_IOS )
  set( INSTALL_DEPENDENCIES_IN_BUILD_DIR OFF )
endif( InterSpec_IOS )

message("")
message("")
message("It may take 10 or 20 minutes to clone into boost and Wt repos...")
message("")


set( BOOST_IOSTREAMS_ENABLE_ZSTD OFF CACHE INTERNAL ""  )
set( BOOST_IOSTREAMS_ENABLE_ZLIB OFF CACHE INTERNAL ""  )
set( BOOST_IOSTREAMS_ENABLE_LZMA OFF CACHE INTERNAL ""  )
set( BOOST_IOSTREAMS_ENABLE_BZIP2 OFF CACHE INTERNAL ""  )
set( Boost_NOWIDE_WERROR OFF CACHE INTERNAL ""  )
#set( BOOST_EXCLUDE_LIBRARIES nowide log CACHE INTERNAL ""  )

set( BUILD_SHARED_LIBS OFF CACHE INTERNAL "Build SHARED libraries" )

# TODO: set BOOST libraries to not build, like wave, beast, wserialization, etc

# CMAKE_CURRENT_SOURCE_DIR could either be project base directory, or a target 
#  directory, so we will search for the patch files to apply
find_file( WT_PATCH_FILE "wt/3.7.1/FetchContent/wt_3.7.1_FetchContent.git.patch"
  PATHS "${CMAKE_CURRENT_SOURCE_DIR}/target/patches"
        "${CMAKE_CURRENT_SOURCE_DIR}/../target/patches"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../target/patches"
  REQUIRED
  NO_DEFAULT_PATH
  NO_CMAKE_FIND_ROOT_PATH
) 

find_file( ZLIB_PATCH_FILE "zlib/1.2.12/FetchContent/zlib_1.2.12.git.patch"  
  PATHS "${CMAKE_CURRENT_SOURCE_DIR}/target/patches"
        "${CMAKE_CURRENT_SOURCE_DIR}/../target/patches"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../target/patches"
  REQUIRED
  NO_DEFAULT_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

find_file( CERES_PATCH_FILE "ceres/2.2.0/ceres.patch"  
  PATHS "${CMAKE_CURRENT_SOURCE_DIR}/target/patches"
        "${CMAKE_CURRENT_SOURCE_DIR}/../target/patches"
        "${CMAKE_CURRENT_SOURCE_DIR}/../../target/patches"
  REQUIRED
  NO_DEFAULT_PATH
  NO_CMAKE_FIND_ROOT_PATH
)

# set( FETCHCONTENT_QUIET FALSE )

# Check for local Boost
set(LOCAL_BOOST_DIR "${CMAKE_BINARY_DIR}/_deps/boost-src")
if(EXISTS "${LOCAL_BOOST_DIR}/CMakeLists.txt")
  message(STATUS "Using local Boost from ${LOCAL_BOOST_DIR}")
  set(BOOST_FETCHCONTENT_SOURCE SOURCE_DIR "${LOCAL_BOOST_DIR}")
  FetchContent_Declare( boost ${BOOST_FETCHCONTENT_SOURCE} )
else()
  set(BOOST_FETCHCONTENT_SOURCE
    #URL https://github.com/boostorg/boost/releases/download/boost-1.85.0/boost-1.85.0-cmake.zip
    #URL_HASH MD5=0f0d58e89e08f8b5f6a01656e19b8535
    #DOWNLOAD_EXTRACT_TIMESTAMP true
    GIT_REPOSITORY https://github.com/boostorg/boost.git
    GIT_TAG        ab7968a0bbcf574a7859240d1d8443f58ed6f6cf # release-1.85.0
    GIT_SHALLOW    ON
  )
  #set(BOOST_ENABLE_LIBRARIES regex thread date_time system filesystem program_options random CACHE STRING "Boost libraries to build")
  #set(Boost_USE_STATIC_LIBS ON CACHE BOOL "Use shared Boost libraries")
  #set(Boost_USE_MULTITHREADED ON CACHE BOOL "Use multithreaded Boost libraries")
  #set(Boost_USE_STATIC_RUNTIME ON CACHE BOOL "Use static runtime for Boost")
  FetchContent_Declare( boost ${BOOST_FETCHCONTENT_SOURCE} )
endif()

FetchContent_MakeAvailable( boost )
FetchContent_GetProperties( boost )

# Check for local Wt
set(LOCAL_WT_DIR "${CMAKE_BINARY_DIR}/_deps/wt-src")
if(EXISTS "${LOCAL_WT_DIR}/CMakeLists.txt")
  message(STATUS "Using local Wt from ${LOCAL_WT_DIR}")
  set(WT_FETCHCONTENT_SOURCE SOURCE_DIR "${LOCAL_WT_DIR}")
else()
  set(WT_FETCHCONTENT_SOURCE
    GIT_REPOSITORY https://github.com/emweb/wt.git
    GIT_TAG        b84925215d2b45879cf20c0cb340c4e7960d0c53 # 3.7.1
    GIT_SHALLOW    ON
    PATCH_COMMAND ${GIT_EXECUTABLE} apply --reverse --check --ignore-space-change --ignore-whitespace ${WT_PATCH_FILE} || ${GIT_EXECUTABLE} apply --reject --ignore-space-change --ignore-whitespace ${WT_PATCH_FILE}
  )
endif()

FetchContent_Declare(
  wt
  ${WT_FETCHCONTENT_SOURCE}
)


if( INSTALL_DEPENDENCIES_IN_BUILD_DIR )
  # We will explicitly populate Wt and boost, rather than using this next line, so this way when the user
  #  builds the install target (e.g., for electron build), it wont install stuff to /usr/local
  FetchContent_MakeAvailable( wt boost )
else( INSTALL_DEPENDENCIES_IN_BUILD_DIR )
  FetchContent_GetProperties(wt)
  if(NOT wt_POPULATED)
    FetchContent_Populate(wt)
    add_subdirectory(${wt_SOURCE_DIR} ${wt_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()

  #FetchContent_GetProperties(boost)
  #if(NOT boost_POPULATED)
  #  FetchContent_Populate(boost)
  #  add_subdirectory(${boost_SOURCE_DIR} ${boost_BINARY_DIR} EXCLUDE_FROM_ALL)
  #endif()

  # Since we arent installing Wt, the InterSpec CMakeLists.txt wont find Wt resources,
  #  so we'll hard code this directory.
  set(WT_RESOURCES_DIRECTORY "${wt_SOURCE_DIR}/resources" CACHE INTERNAL "Set Wt resources dir")

  # If we run the Cmake "install" target (e.g., for Electron build), we'll get an error because 
  # ${boost_BINARY_DIR}/cmake_install.cmake doesnt exist, so we'll just put a dummy 
  # there, and things seem to work
  #   However, the proper fix would be something like:
  #     https://stackoverflow.com/questions/65527126/disable-install-for-fetchcontent#answer-65529161
  #   Or maybe even better yet, we should modify the patching for Wt and zlib to remove the instal(...) commands
  write_file( "${boost_BINARY_DIR}/cmake_install.cmake" "# This is a dummy patch file the InterSpec build put in" )
endif( INSTALL_DEPENDENCIES_IN_BUILD_DIR )



# macOS, iOS, and Android all have zlib already - we dont need to build it
if( NOT APPLE AND NOT ANDROID )
  set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "Build SHARED libraries")

  FetchContent_Declare(
    zlib
    GIT_REPOSITORY https://github.com/madler/zlib.git
    GIT_TAG        21767c654d31d2dccdde4330529775c6c5fd5389 # 1.2.12
    GIT_SHALLOW    ON

    PATCH_COMMAND ${GIT_EXECUTABLE} apply --reverse --check --ignore-space-change --ignore-whitespace ${ZLIB_PATCH_FILE} || ${GIT_EXECUTABLE} apply --reject --ignore-space-change --ignore-whitespace ${ZLIB_PATCH_FILE}
    #PATCH_COMMAND patch -p1 < ${ZLIB_PATCH_FILE} || true
  )

  FetchContent_GetProperties(zlib)
  if(NOT zlib_POPULATED)
    FetchContent_Populate(zlib)
    add_subdirectory(${zlib_SOURCE_DIR} ${zlib_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()

  FetchContent_MakeAvailable( zlib )
   
  if( ZLIB::zlib ) 
    # Hmmm, not sure if we ever actually get here
    add_library( ZLIB::ZLIB ALIAS ZLIB::zlib )
  elseif( NOT ZLIB::ZLIB )
    # On Linux, we seem to always get here, so we'll just hack things.  
    #  I must not understand the export(...) cmake stuff of the sub-project - oh well for now.
    add_library( ZLIB::ZLIB ALIAS zlib )
    set( ZLIB_LIBRARY zlib CACHE INTERNAL "set zlib static lib" )
    set( ZLIB_INCLUDE_DIR "${zlib_BINARY_DIR}" CACHE INTERNAL "set zlib include dir" )
  endif()
endif( NOT APPLE AND NOT ANDROID )


  
# Set some eigen options
set( EIGEN_BUILD_SHARED_LIBS OFF CACHE INTERNAL "" )
set( EIGEN_BUILD_DOC OFF CACHE INTERNAL "" )
set( EIGEN_BUILD_TESTING OFF CACHE INTERNAL "" )
set( EIGEN_BUILD_PKGCONFIG OFF )
set( EIGEN_BUILD_CMAKE_PACKAGE ON CACHE INTERNAL "" )
#set( Eigen3_DIR "${FETCHCONTENT_BASE_DIR}/..." CACHE INTERNAL "" )
# Check for local Eigen
set(LOCAL_EIGEN_DIR "${CMAKE_BINARY_DIR}/_deps/eigen-src")
if(EXISTS "${LOCAL_EIGEN_DIR}/CMakeLists.txt")
  message(STATUS "Using local Eigen from ${LOCAL_EIGEN_DIR}")
  set(EIGEN_FETCHCONTENT_SOURCE SOURCE_DIR "${LOCAL_EIGEN_DIR}")
else()
  set(EIGEN_FETCHCONTENT_SOURCE
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        2e76277bd049f7bec36b0f908c69734a42c5234f  # Get trunk version as of 20250114 to pickup some compile issues for c++20
    UPDATE_DISCONNECTED TRUE # Prevent automatic updates if the tag is already present
  )
endif(EXISTS "${LOCAL_EIGEN_DIR}/CMakeLists.txt")

FetchContent_Declare( eigen ${EIGEN_FETCHCONTENT_SOURCE} )
FetchContent_MakeAvailable( eigen )
#FetchContent_GetProperties( eigen )
#if(NOT eigen_POPULATED)
#  FetchContent_Populate(eigen)
#  add_subdirectory(${eigen_SOURCE_DIR} ${eigen_BINARY_DIR} EXCLUDE_FROM_ALL)
#endif()
  
  # For Android and iOS, we need to force the path information for Eigen, for some reason.
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${eigen_BINARY_DIR}" CACHE STRING "Modules for CMake" FORCE )
set( Eigen3_DIR "${eigen_BINARY_DIR}/share/eigen3/cmake/" CACHE PATH "Path to Eigen3Config.cmake" )  #${eigen_SOURCE_DIR}

  # Set some Google Ceres options; both to allow compiling without further dependencies, and
  #  avoid some things we dont need
set( MINIGLOG ON CACHE INTERNAL "" )
set( GFLAGS OFF CACHE INTERNAL "" )
set( ACCELERATESPARSE OFF CACHE INTERNAL "" )
set( USE_CUDA OFF CACHE INTERNAL "" )
set( EXPORT_BUILD_DIR OFF CACHE INTERNAL "" )
set( BUILD_TESTING OFF CACHE INTERNAL "" )
set( BUILD_EXAMPLES OFF CACHE INTERNAL "" )
set( PROVIDE_UNINSTALL_TARGET OFF CACHE INTERNAL "" )
# Check for local Ceres Solver
set(LOCAL_CERES_DIR "${CMAKE_BINARY_DIR}/_deps/ceres-solver-src")
if(EXISTS "${LOCAL_CERES_DIR}/CMakeLists.txt")
  message(STATUS "Using local Ceres Solver from ${LOCAL_CERES_DIR}")
  set(CERES_FETCHCONTENT_SOURCE SOURCE_DIR "${LOCAL_CERES_DIR}" EXCLUDE_FROM_ALL)
else()
  set(CERES_FETCHCONTENT_SOURCE
      GIT_REPOSITORY https://github.com/ceres-solver/ceres-solver.git
      GIT_TAG         85331393dc0dff09f6fb9903ab0c4bfa3e134b01 # release-2.2.0
      UPDATE_DISCONNECTED TRUE
      GIT_SHALLOW    TRUE
      PATCH_COMMAND ${GIT_EXECUTABLE} apply --reverse --check --ignore-space-change --ignore-whitespace ${CERES_PATCH_FILE} || ${GIT_EXECUTABLE} apply --reject --ignore-space-change --ignore-whitespace ${CERES_PATCH_FILE}
      EXCLUDE_FROM_ALL
  )
endif()

FetchContent_Declare( ceres-solver ${CERES_FETCHCONTENT_SOURCE} )
FetchContent_MakeAvailable( ceres-solver )
#FetchContent_GetProperties( ceres-solver )
#if(NOT ceres-solver_POPULATED)
#  FetchContent_Populate(ceres-solver)
#  add_subdirectory(${ceres-solver_SOURCE_DIR} ${ceres_solver_BINARY_DIR} EXCLUDE_FROM_ALL)
#endif()

if( BUILD_AS_WX_WIDGETS_APP )
  set(wxBUILD_SHARED OFF CACHE BOOL "")

  if(WIN32)
    set(wxUSE_WEBVIEW_EDGE ON CACHE BOOL "")
    set(wxUSE_WEBVIEW_EDGE_STATIC ON CACHE BOOL "")
    # TODO: static runtime should be set off of wxInterSpec_USE_MSVC_MultiThreadDLL
    set(wxBUILD_USE_STATIC_RUNTIME ON CACHE BOOL "")
  endif(WIN32)

  # Think by default sub-modules are also retrieved - but not tested yet
  FetchContent_Declare( wxWidgets
    GIT_REPOSITORY https://github.com/wxWidgets/wxWidgets.git
    GIT_TAG        97e99707c5d2271a70cb686720b48dbf34ced496 # release-3.2.1
    GIT_SHALLOW    ON
    # GIT_PROGRESS TRUE # Also need FETCHCONTENT_QUIET 
  )

  FetchContent_GetProperties(wxWidgets)
  if(NOT wxWidgets_POPULATED)
    FetchContent_Populate(wxWidgets)
    add_subdirectory(${wxWidgets_SOURCE_DIR} ${wxWidgets_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif(NOT wxWidgets_POPULATED)

  # write_file( "${wxWidgets_BINARY_DIR}/cmake_install.cmake" "# This is a dummy patch file the InterSpec build put in" )
endif( BUILD_AS_WX_WIDGETS_APP )

message("Building Wt and boost may also take up to a couple hours...")