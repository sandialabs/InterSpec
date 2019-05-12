project(InterSpec)

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif( COMMAND cmake_policy )


set( CMAKE_COLOR_MAKEFILE ON CACHE BOOL "Color" FORCE )

set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

IF(WIN32)
  set( MSVC_RUNTIME "static" )
	set( TRY_TO_STATIC_LINK "ON" )
  include( cmake/ConfigureMsvc.txt )
	configure_msvc_runtime() 

  set( BOOST_PREFIX "C:/install/msvc2017/x64/boost_1_65_1" )
  set( ZLIB_PREFIX "C:/install/msvc2017/x64/zlib" )
  list( APPEND PATH_SUFFIXES "C:/install/msvc2017/x64/wt-3.3.4" )
  
  set( EXECUTABLE_OUTPUT_PATH . CACHE PATH "Path to executables" FORCE )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj -D_SCL_SECURE_NO_WARNINGS /MP /wd4996 /wd4267 /DWINVER=0x0601 /D_WIN32_WINNT=0x0601" )
    #0x0601==Win7, 0x0501==WinXP
ELSE(WIN32)
    set( EXECUTABLE_OUTPUT_PATH bin CACHE PATH "Path to executables" FORCE )
#    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=512" )
ENDIF(WIN32)

SET( InterSpec_VERSION Development CACHE STRING "InterSpec Version" )
SET( MYSQL_DATABASE_TO_USE "dev" CACHE STRING "Which mysql database to use (if applicable): dev/qc/prod" )

option( BUILD_AS_OSX_APP "Build a OSX native app" OFF )
option( BUILD_AS_ELECTRON_APP "Configures the executable to be ran as a Electron app" OFF )
option( BUILD_FOR_WEB_DEPLOYMENT "Selects options appropriate for the web in the future" OFF )
option( BUILD_AS_LOCAL_SERVER "Build for use as a server for local computer (e.g. run from the command line and connect to with a browser on LOCALHOST)" ON )
option( BUILD_AS_UNIT_TEST_SUITE "Builds unit tests as well the analysis tests (aka end-to-end tests)" OFF )
option( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE "Compiles so executable only does offline testing of the user test states in the database (e.g. end-to-end testing)" OFF )
option( BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT "Build executable for use while testing new code, not running InterSpec" OFF )

option( INCLUDE_ANALYSIS_TEST_SUITE "Allow whether user can save and load test spectra" on )
set( TEST_SUITE_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/testing" CACHE STRING "Path to directory that contains the \"analysis_tests\" directory for saving N42 test states.  Leave empty for CWD." )

option( PERFORM_DEVELOPER_CHECKS "Performs additional computationally expensive tests during execution" off )
option( USE_OSX_NATIVE_MENU "Mirrors WMenus with OSX native NSMenu (not fully implemented)" off )
option( USE_ELECTRON_NATIVE_MENU "Mirrors WMenus with Electrons native Menu (not fully implemented)" on )
option( USE_HIGH_BANDWIDTH_INTERACTIONS "Allow more interactiveness with the charts" ON )
option( SUPPORT_ZIPPED_SPECTRUM_FILES "(Work in progress) Allows opening zipped spectrum files" ON )
option( ENABLE_RESOURCE_UPDATES "Enables user uploads of resource (nuc decay info, DRFs, etc) upgrades" ON )
option( TRY_TO_STATIC_LINK "Try to link to static libs whenever possible" OFF )
option( DRAW_GAMMA_LINES_LOG_AND_LIN "Have gamma lines indicators be drawn in log mode when y-axis is log" OFF )
option( ALLOW_URL_TO_FILESYSTEM_MAP "Allows a file ID to be specified in the URL to open a file from the filesystem" ON )
option( USE_DB_TO_STORE_SPECTRA "Use the database to store spectra" ON )
option( RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE "Render photopeaks server-side rather than client side" ON )
option( USE_SPECRUM_FILE_QUERY_WIDGET "Enables compilation of the spectrum file query tool (not for web deployment)" OFF )
option( USE_SEARCH_MODE_3D_CHART "Enable 3D chart or searchmode/passthrough data.  Requires Wt >=3.3.4 compiled with WebGL/OpenGL support" OFF )
option( USE_TERMINAL_WIDGET "Enables use of terminal widget" OFF )
set( SpecUtils_ENABLE_D3_CHART ON CACHE BOOL "Enables use of exporting D3 HTML files" )
option( USE_SPECTRUM_CHART_D3 "(experimental) Replaces Wt SpectrumChart with D3.js chart" OFF )
option( USE_SIMPLE_NUCLIDE_ASSIST "Enable using of Mike E. gamma rays of interest to help with nuclide ID" OFF )

option( USE_SQLITE3_DB "Use SQLITE3 database to store user preferences and spectra" ON )
option( USE_MYSQL_DB "Use MySQL database to store user preferences and spectra" OFF )
option( USE_GOOGLE_MAP "Use google maps widget" ON )
option( SpecUtils_REBIN_FILES_TO_SINGLE_BINNING "Rebins all spectra in spectrum file to the same energy calibration" ON )
option( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK "Include support for downloading images of the displayed page" OFF )

set( MAX_SPECTRUM_MEMMORY_SIZE_MB 256 CACHE STRING "Amount of memory to allow spectra to take up before trying to offload them onto disk when not in use" )

set( GOOGLE_MAPS_KEY "" CACHE STRING "Google maps api key." )

set( GUI_TYPE "" )

IF( WIN32 OR BUILD_AS_OSX_APP OR BUILD_AS_ELECTRON_APP )
    set( TRY_TO_STATIC_LINK "ON" )
ENDIF( WIN32 OR BUILD_AS_OSX_APP OR BUILD_AS_ELECTRON_APP )


IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(POLICY CMP0048)  
 	  CMAKE_POLICY(SET CMP0048 OLD)
ENDIF(POLICY CMP0048)

if(BUILD_AS_UNIT_TEST_SUITE)
    set(BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE ON)
endif(BUILD_AS_UNIT_TEST_SUITE)

if(BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE)
    set(INCLUDE_ANALYSIS_TEST_SUITE ON)
endif(BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE)

if(INCLUDE_ANALYSIS_TEST_SUITE)
    set(USE_DB_TO_STORE_SPECTRA ON)
endif(INCLUDE_ANALYSIS_TEST_SUITE)

if(BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)
    set(PERFORM_DEVELOPER_CHECKS ON)
endif(BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)


if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
   #get rid of about a million warnings caused by boost foreach
   add_definitions ("-Wno-unused-local-typedefs")
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   #get rid of a million warnings by Wt/Dbo/ptr_impl.h in VS2015
   add_definitions( "/wd4297" )
   add_definitions( "/wd4244" )
   #set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/ignore:4099")
endif()


#-----------------------------------------------------------------------------

if( BUILD_AS_OSX_APP )
    SET(EXECUTABLE_NAME InterSpec)
else( BUILD_AS_OSX_APP )
    SET(LIBRARYNAME LibInterSpec)
    IF(WIN32)
        SET(EXECUTABLE_NAME InterSpec)
    else(WIN32)
        SET(EXECUTABLE_NAME InterSpec.exe)
    endif(WIN32)
    if(ANDROID)
      SET( LIBRARYNAME "InterSpec" )
    endif()
endif( BUILD_AS_OSX_APP )

if(ANDROID OR IOS OR BUILD_AS_UNIT_TEST_SUITE )
  SET( LINK_TO_OBJECT ${LIBRARYNAME})
else(ANDROID OR IOS OR BUILD_AS_UNIT_TEST_SUITE )
  SET( LINK_TO_OBJECT ${EXECUTABLE_NAME})
endif(ANDROID OR IOS OR BUILD_AS_UNIT_TEST_SUITE )


IF(ANDROID)
  set(BUILD_AS_LOCAL_SERVER OFF)
  set(ALLOW_URL_TO_FILESYSTEM_MAP ON)
  set(TRY_TO_STATIC_LINK ON)
  set(MAX_SPECTRUM_MEMMORY_SIZE_MB 32)
  set(USE_DB_TO_STORE_SPECTRA ON)
  set(USE_HIGH_BANDWIDTH_INTERACTION ON)
  set(USE_SPECRUM_FILE_QUERY_WIDGET OFF)
  set(USE_TERMINAL_WIDGET OFF)
  set(SpecUtils_ENABLE_D3_CHART OFF)
  set(USE_SPECTRUM_CHART_D3 OFF)
  set(USE_SEARCH_MODE_3D_CHART ON)
  set(INCLUDE_ANALYSIS_TEST_SUITE OFF)
  set(SHARED_LIBS OFF)
  add_definitions(-DANDROID -DWT_NO_STD_WSTRING=ON -DWT_NO_STD_LOCALE=ON )


#This is just painful.  I cant get the find_* macros to work with the CMake
# included as part of Android Studio... After wasting what seems like an entire
# day, its time for the brute force solution
  if( NOT BOOST_WT_PREFIX_BASE )
   message( fatal "You must define BOOST_WT_PREFIX_BASE when building Android build" )
  endif( NOT BOOST_WT_PREFIX_BASE )
  set( USRPREFIX "${BOOST_WT_PREFIX_BASE}/${ANDROID_ABI}" )
  set( BOOST_LIB_START "${USRPREFIX}/lib/libboost_" )
  set( BOOST_LIB_END "-clang-darwin-mt-1_65_1.a" )

  list( APPEND CMAKE_FIND_ROOT_PATH "${USRPREFIX}" )

  set( BOOST_PREFIX "${USRPREFIX}" )
  set( Boost_LIBRARY_DIR "${USRPREFIX}/lib" )
  set( Boost_INCLUDE_DIR "${USRPREFIX}/include" )
  set( Boost_ATOMIC_LIBRARY_RELEASE ${BOOST_LIB_START}atomic${BOOST_LIB_END} )
  set( Boost_CHRONO_LIBRARY_RELEASE ${BOOST_LIB_START}chrono${BOOST_LIB_END} )
  set( Boost_DATE_TIME_LIBRARY_RELEASE ${BOOST_LIB_START}date_time${BOOST_LIB_END} )
  set( Boost_FILESYSTEM_LIBRARY_RELEASE ${BOOST_LIB_START}filesystem${BOOST_LIB_END} )
  set( Boost_IOSTREAMS_LIBRARY_RELEASE ${BOOST_LIB_START}iostreams${BOOST_LIB_END} )
  set( Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE ${BOOST_LIB_START}program_options${BOOST_LIB_END} )
  set( Boost_RANDOM_LIBRARY_RELEASE ${BOOST_LIB_START}random${BOOST_LIB_END} )
  set( Boost_REGEX_LIBRARY_RELEASE ${BOOST_LIB_START}regex${BOOST_LIB_END} )
  set( Boost_SYSTEM_LIBRARY_RELEASE ${BOOST_LIB_START}system${BOOST_LIB_END} )
  set( Boost_THREAD_LIBRARY_RELEASE ${BOOST_LIB_START}thread${BOOST_LIB_END} )
  set( Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE ${BOOST_LIB_START}unit_test_framework${BOOST_LIB_END} )
  set( Boost_SIGNALS_LIBRARY_RELEASE ${BOOST_LIB_START}signals${BOOST_LIB_END} )
  set( Boost_SERIALIZATION_LIBRARY_RELEASE ${BOOST_LIB_START}serialization${BOOST_LIB_END} )
  set( Wt_INCLUDE_DIR=${USRPREFIX}/include )
  set( Wt_LIBRARY=${USRPREFIX}/lib/libwt.a )
  set( Wt_TEST_LIBRARY=${USRPREFIX}/lib/libwttest.a )
  set( Wt_HTTP_LIBRARY=${USRPREFIX}/lib/libwthttp.a )
  set( Wt_DBO_LIBRARY=${USRPREFIX}/lib/libwtdbo.a )
  set( Wt_DBOSQLITE3_LIBRARY=${USRPREFIX}/lib/libwtdbosqlite3.a )
ENDIF(ANDROID)


IF(IOS)
  set( CMAKE_XCODE_ATTRIBUTE_IPHONEOS_DEPLOYMENT_TARGET "7.0" CACHE STRING "Set the minimum deployment target value." FORCE)

  #For iOS we have to use libtool instead of ar so that the platform (armv7, armv7s, arm64, etc) are kept track of in the archive
  SET(CMAKE_AR "libtool" )

  #CMake adds a cq after the archive command by default, we need to change this to "-static -o"
  SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> -static -o <TARGET> <LINK_FLAGS> <OBJECTS>")
  SET(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> -static -o <TARGET> <LINK_FLAGS> <OBJECTS>")

  set(BUILD_AS_LOCAL_SERVER OFF)
  set(ALLOW_URL_TO_FILESYSTEM_MAP ON)
  set(INCLUDE_ANALYSIS_TEST_SUITE OFF)
  set(PERFORM_DEVELOPER_CHECKS OFF)
  set(TRY_TO_STATIC_LINK ON)
  set(USE_SPECRUM_FILE_QUERY_WIDGET OFF)
  set(MAX_SPECTRUM_MEMMORY_SIZE_MB "32" )
  set(FRAMEWORKDIR "${CMAKE_CURRENT_SOURCE_DIR}/target/ios")
  set(CMAKE_SYSTEM_FRAMEWORK_PATH "${FRAMEWORKDIR}")
  set(USE_BOOST_FRAMEWORK OFF)
  set(DEFAULT_WT_BOOST_DISCOVERY TRUE)
  set(CMAKE_FRAMEWORK_PATH "${FRAMEWORKDIR}")
  set(BOOST_FRAMEWORK_PATH "${CMAKE_FRAMEWORK_PATH}/framework/boost.framework")
  set(BOOST_PREFIX "${FRAMEWORKDIR}/prefix")
  set(Boost_DIR "${FRAMEWORKDIR}/prefix" )
  set(Boost_INCLUDE_DIR "${Boost_DIR}/include" )
  set(Boost_LIBRARY_DIRS "${Boost_DIR}/lib" )
  set(BOOST_LIBRARYDIR "${Boost_DIR}/lib" )
  set(BOOST_ROOT "${Boost_DIR}" )
  SET(CMAKE_MODULE_PATH ${Boost_DIR} ${CMAKE_MODULE_PATH} )
  set(SHARED_LIBS OFF)
  set(CMAKE_FIND_ROOT_PATH ${CMAKE_FIND_ROOT_PATH} ${Boost_DIR})
  set(SUPPORT_ZIPPED_SPECTRUM_FILES ON)
  set(USE_TERMINAL_WIDGET OFF)
  set(SpecUtils_ENABLE_D3_CHART OFF)
  set(USE_SPECTRUM_CHART_D3 OFF)

  #if you change this IPHONEOS_DEPLOYMENT_TARGET also change the variable in
  # iosbuild.rb everywhere, and in the Xcode project.
  set(ENV{IPHONEOS_DEPLOYMENT_TARGET} 10.2)

  #We know where Wt is for iOS, so just force it
  IF ( NOT Wt_INCLUDE_DIR )
    SET (Wt_DIR "${CMAKE_CURRENT_SOURCE_DIR}/target/ios/prefix")
    SET (Wt_INCLUDE_DIR "${Wt_DIR}/include")
  ENDIF ( NOT Wt_INCLUDE_DIR )
ENDIF(IOS)

IF( BUILD_AS_OSX_APP )
  set( ALLOW_URL_TO_FILESYSTEM_MAP ON )
  set( USE_OSX_NATIVE_MENU ON )
  set( USE_SPECRUM_FILE_QUERY_WIDGET ON )
  set( USE_TERMINAL_WIDGET OFF )
  set( SpecUtils_ENABLE_D3_CHART ON )
  set( USE_SPECTRUM_CHART_D3 OFF )
ENDIF( BUILD_AS_OSX_APP )


IF( BUILD_AS_ELECTRON_APP )
  set( BUILD_AS_LOCAL_SERVER OFF )
  set( ALLOW_URL_TO_FILESYSTEM_MAP ON )
  set( USE_OSX_NATIVE_MENU OFF )
  set( USE_SPECRUM_FILE_QUERY_WIDGET ON )
  set( USE_TERMINAL_WIDGET OFF )
ENDIF( BUILD_AS_ELECTRON_APP )


IF( BUILD_AS_UNIT_TEST_SUITE )
    set( BUILD_AS_LOCAL_SERVER OFF )
    set( USE_ENERGY_VS_TIME_PLOTS OFF )
    set( USE_DB_TO_STORE_SPECTRA ON )
    set( ALLOW_URL_TO_FILESYSTEM_MAP ON )
    set( INCLUDE_ANALYSIS_TEST_SUITE ON )
    set( TRY_TO_STATIC_LINK ON )
    set( PERFORM_DEVELOPER_CHECKS ON )
    set( USE_SPECRUM_FILE_QUERY_WIDGET ON )
    set( USE_TERMINAL_WIDGET ON )
    set( SpecUtils_ENABLE_D3_CHART ON )
    set( USE_SPECTRUM_CHART_D3 OFF )
ENDIF( BUILD_AS_UNIT_TEST_SUITE )

IF( BUILD_FOR_WEB_DEPLOYMENT )
    set( TRY_TO_STATIC_LINK ON )
    set( CMAKE_BUILD_TYPE:STRING "Release" )
    set( ALLOW_URL_TO_FILESYSTEM_MAP OFF )
    set( USE_DB_TO_STORE_SPECTRA ON )
#    set( USE_SQLITE3_DB ON )
#    set( USE_MYSQL_DB OFF )
    set( DATABASE_PASSWORD_FILE "/path/to/passwords/databases.xml" )
    set( USE_GOOGLE_MAP OFF )
    set( INCLUDE_ANALYSIS_TEST_SUITE OFF )
    set( PERFORM_DEVELOPER_CHECKS OFF )
    set( SUPPORT_ZIPPED_SPECTRUM_FILES ON )
    set( USE_SPECRUM_FILE_QUERY_WIDGET OFF )
    set( USE_TERMINAL_WIDGET OFF )
    set( SpecUtils_ENABLE_D3_CHART ON )
    set( USE_SPECTRUM_CHART_D3 OFF )
    set( ENABLE_RESOURCE_UPDATES OFF )
ENDIF( BUILD_FOR_WEB_DEPLOYMENT )


#check to make sure exactly one target is specified
set( NUM_TARGETS_SPECIFIED 0 )

#if( BUILD_AS_OSX_APP )
#  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
#endif()

if( BUILD_AS_ELECTRON_APP )
  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if( BUILD_FOR_WEB_DEPLOYMENT )
  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if( BUILD_AS_LOCAL_SERVER )
  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE )
  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if( BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )
  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if( ANDROID OR IOS )
  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if( NOT  ${NUM_TARGETS_SPECIFIED} EQUAL "1" )
  message( FATAL "You must specify exactly one target, you specified ${NUM_TARGETS_SPECIFIED}" )
endif()


#try to link to static libraries whenever possible
#  doing this to boost libraries only adds about 0.6 megabytes
if(TRY_TO_STATIC_LINK)
    set(Boost_USE_STATIC_LIBS ON)
    IF(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    ELSE(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    ENDIF(WIN32)
else(TRY_TO_STATIC_LINK)
    set(Boost_USE_STATIC_LIBS OFF)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif(TRY_TO_STATIC_LINK)


INCLUDE(cmake/FindWt.cmake)
INCLUDE(cmake/WtFindZlib.txt)  #ToDo: switch to using CMake privided FindZLIB.cmake, e.g. use find_package( ZLIB REQUIRED )

SET( Boost_ADDITIONAL_VERSIONS "1.55" "1.57" "1.65" "1.65_1" )
#find_package( Boost 1.41.0 COMPONENTS thread date_time system filesystem program_options regex random chrono atomic )

INCLUDE(cmake/WtFindBoost-cmake.txt)
IF(NOT BOOST_WT_FOUND)
  SET(ERR
      "Could not find a boost installation in " ${BOOST_PREFIX} ".\n\n"
      "There are two methods in Wt to find boost:")
  MESSAGE(FATAL_ERROR ${ERR})
ENDIF(NOT BOOST_WT_FOUND)


set( sources
    src/DecayDataBaseServer.cpp
    src/IsotopeSelectionAids.cpp
    src/IsotopeId.cpp
    src/MaterialDB.cpp
    src/PhysicalUnits.cpp
    src/MassAttenuationTool.cpp
    src/DetectorPeakResponse.cpp
    src/GammaInteractionCalc.cpp
    src/ShieldingSourceDisplay.cpp
    src/GammaXsGui.cpp
    src/ReferenceLineInfo.cpp
    src/ReferencePhotopeakDisplay.cpp
    src/IsotopeNameFilterModel.cpp
    src/ReactionGamma.cpp
    src/IsotopeSearchByEnergy.cpp
    src/IsotopeSearchByEnergyModel.cpp
    src/SpectrumChart.cpp
    src/SpectrumDataModel.cpp
    src/SpectrumDisplayDiv.cpp
    src/HelpSystem.cpp
    src/InterSpec.cpp
    src/PopupDiv.cpp
    src/CanvasForDragging.cpp
    src/SpecMeas.cpp
    src/PeakFit.cpp
    src/PeakDef.cpp
    src/SpectraFileModel.cpp
    src/AuxWindow.cpp
    src/PeakFitChi2Fcn.cpp
    src/PeakModel.cpp
    src/PeakInfoDisplay.cpp
    src/WarningWidget.cpp
    src/SpecMeasManager.cpp
    src/Recalibrator.cpp
    src/InterSpecUser.cpp
    src/CompactFileManager.cpp
    src/DetectorEdit.cpp
    src/OneOverR2Calc.cpp
    src/ActivityConverter.cpp
    src/DecayWindow.cpp
    src/DecayActivityDiv.cpp
    src/DecayChainChart.cpp
    src/DecaySelectNuclideDiv.cpp
    src/SpecFileSummary.cpp
    src/GammaCountDialog.cpp
    src/PeakEdit.cpp
    src/UseInfoWindow.cpp
    src/LicenseAndDisclaimersWindow.cpp
    src/LocalTimeDelegate.cpp
    src/RowStretchTreeView.cpp
    src/DataBaseVersionUpgrade.cpp
    src/DataBaseUtils.cpp
    src/Integrate.cpp
    src/GadrasSpecFunc.cpp
    src/DoseCalc.cpp
    src/DoseCalcWidget.cpp
    src/GadrasGamFileParser.cpp
    src/ColorSelect.cpp
    src/ColorTheme.cpp
    src/ColorThemeWidget.cpp
    src/ColorThemeWindow.cpp
    src/PeakSearchGuiUtils.cpp
    src/MakeDrf.cpp
    src/MakeDrfSrcDef.cpp
    src/MakeDrfChart.cpp
    src/MakeDrfFit.cpp
    js/CanvasForDragging.js
    js/SpectrumChart.js
    js/InterSpec.js
)


set( headers
    InterSpec/InterSpec_config.h.in
    InterSpec/DecayDataBaseServer.h
    InterSpec/IsotopeSelectionAids.h
    InterSpec/IsotopeId.h
    InterSpec/MaterialDB.h
    InterSpec/PhysicalUnits.h
    InterSpec/MassAttenuationTool.h
    InterSpec/DetectorPeakResponse.h
    InterSpec/GammaInteractionCalc.h
    InterSpec/ShieldingSourceDisplay.h
    InterSpec/GammaXsGui.h
    InterSpec/ReferenceLineInfo.h
    InterSpec/ReferencePhotopeakDisplay.h
    InterSpec/IsotopeNameFilterModel.h
    InterSpec/ReactionGamma.h
    InterSpec/IsotopeSearchByEnergy.h
    InterSpec/IsotopeSearchByEnergyModel.h
    InterSpec/SpectrumChart.h
    InterSpec/SpectrumDataModel.h
    InterSpec/SpectrumDisplayDiv.h
    InterSpec/HelpSystem.h
    InterSpec/InterSpec.h
    InterSpec/PopupDiv.h
    InterSpec/CanvasForDragging.h
    InterSpec/SpecMeas.h
    InterSpec/PeakFit.h
    InterSpec/PeakDef.h
    InterSpec/SpectraFileModel.h
    InterSpec/AuxWindow.h
    InterSpec/PeakFitChi2Fcn.h
    InterSpec/PeakModel.h
    InterSpec/PeakInfoDisplay.h
    InterSpec/WarningWidget.h
    InterSpec/SpecMeasManager.h
    InterSpec/Recalibrator.h
    InterSpec/InterSpecUser.h
    InterSpec/CompactFileManager.h
    InterSpec/DetectorEdit.h
    InterSpec/OneOverR2Calc.h
    InterSpec/ActivityConverter.h
    InterSpec/DecayWindow.h
    InterSpec/DecayActivityDiv.h
    InterSpec/DecayChainChart.h
    InterSpec/DecaySelectNuclideDiv.h
    InterSpec/SpecFileSummary.h
    InterSpec/GammaCountDialog.h
    InterSpec/PeakEdit.h
    InterSpec/UseInfoWindow.h
    InterSpec/LicenseAndDisclaimersWindow.h
    InterSpec/LocalTimeDelegate.h
    InterSpec/RowStretchTreeView.h
    InterSpec/DataBaseVersionUpgrade.h
    InterSpec/DataBaseUtils.h
    InterSpec/Integrate.h
    InterSpec/GadrasSpecFunc.h
    InterSpec/DoseCalc.h
    InterSpec/DoseCalcWidget.h
    InterSpec/GadrasGamFileParser.h
    InterSpec/ColorSelect.h
    InterSpec/ColorTheme.h
    InterSpec/ColorThemeWidget.h
    InterSpec/ColorThemeWindow.h
    InterSpec/PeakSearchGuiUtils.h
    InterSpec/MakeDrf.h
    InterSpec/MakeDrfSrcDef.h
    InterSpec/MakeDrfChart.h
    InterSpec/MakeDrfFit.h
)

if( USE_DB_TO_STORE_SPECTRA )
  set( headers ${headers} InterSpec/DbStateBrowser.h InterSpec/DbFileBrowser.h )
  set( sources ${sources} src/DbStateBrowser.cpp src/DbFileBrowser.cpp )
endif( USE_DB_TO_STORE_SPECTRA )

if( ALLOW_URL_TO_FILESYSTEM_MAP )
    set( headers ${headers} InterSpec/DbToFilesystemLink.h  )
    set( sources ${sources} src/DbToFilesystemLink.cpp  )
endif( ALLOW_URL_TO_FILESYSTEM_MAP )

if( NOT ANDROID AND NOT IOS )
    set( sources ${sources} src/FileDragUploadResource.cpp )
    set( headers ${headers} InterSpec/FileDragUploadResource.h )
endif( NOT ANDROID AND NOT IOS )

if( BUILD_AS_ELECTRON_APP )
  set( sources ${sources} src/InterSpecServer.cpp target/electron/ElectronUtils.cpp )
  set( headers ${headers} InterSpec/InterSpecServer.h target/electron/ElectronUtils.h )
endif( BUILD_AS_ELECTRON_APP )


set( sources src/InterSpecApp.cpp ${sources} )
set( headers InterSpec/InterSpecApp.h ${headers} )


if( INCLUDE_ANALYSIS_TEST_SUITE )
    set( sources ${sources} src/SpectrumViewerTester.cpp )
    set( headers ${headers} InterSpec/SpectrumViewerTester.h )
endif( INCLUDE_ANALYSIS_TEST_SUITE )

if( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )
  set( SRCS ${SRCS} src/ChartToImageResource.cpp )
  set( HEADRS ${HEADRS} InterSpec/ChartToImageResource.h )
endif( DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK )


if(BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)
    set( sources ${sources} testing/developcode.cpp src/DbUserState.cpp )
    set( headers ${headers} testing/developcode.h InterSpec/DbUserState.h )
endif(BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)


if( USE_OSX_NATIVE_MENU )
    set( sources ${sources} target/osx/NativeMenu.mm )
    set( headers ${headers} target/osx/NativeMenu.h )
endif( USE_OSX_NATIVE_MENU )

if( IOS )
    set( sources ${sources} target/ios/InterSpec/FileHandling.mm )
    set( headers ${headers} target/ios/InterSpec/FileHandling.h )
endif( IOS )

if( USE_GOOGLE_MAP )
    set( sources ${sources} src/GoogleMap.cpp )
    set( headers ${headers} InterSpec/GoogleMap.h )
endif( USE_GOOGLE_MAP )


if( USE_SEARCH_MODE_3D_CHART )
    set( sources ${sources} src/SearchMode3DChart.cpp src/SearchMode3DDataModel.cpp )
    set( headers ${headers} InterSpec/SearchMode3DChart.h InterSpec/SearchMode3DDataModel.h )
endif( USE_SEARCH_MODE_3D_CHART )


if( USE_TERMINAL_WIDGET )
    set( sources ${sources} src/TerminalModel.cpp    src/TerminalWidget.cpp   js/TerminalWidget.js )
    set( headers ${headers} InterSpec/TerminalModel.h   InterSpec/TerminalWidget.h )
endif( USE_TERMINAL_WIDGET )


if( SUPPORT_ZIPPED_SPECTRUM_FILES )
    set( sources ${sources} src/ZipArchive.cpp )
    set( headers ${headers} InterSpec/ZipArchive.h )
endif( SUPPORT_ZIPPED_SPECTRUM_FILES )


if( ENABLE_RESOURCE_UPDATES )
    set( sources ${sources} src/ResourceUpdate.cpp src/ResourceUpdateTool.cpp )
    set( headers ${headers} InterSpec/ResourceUpdate.h InterSpec/ResourceUpdateTool.h )

  if( BUILD_FOR_WEB_DEPLOYMENT )
    MESSAGE( FATAL_ERROR "ENABLE_RESOURCE_UPDATES must not be enabled if you are building for web deployment" )
  endif( BUILD_FOR_WEB_DEPLOYMENT )
endif( ENABLE_RESOURCE_UPDATES )

if( USE_SIMPLE_NUCLIDE_ASSIST )
    set( sources ${sources} src/SimpleNuclideAssist.cpp )
    set( headers ${headers} InterSpec/SimpleNuclideAssist.h )
endif( USE_SIMPLE_NUCLIDE_ASSIST )


IF( USE_SPECRUM_FILE_QUERY_WIDGET )
  IF( BUILD_FOR_WEB_DEPLOYMENT )
    message ( FATAL_ERROR "The spectrum file query widget can not be enabled for web deployements." )
  ENDIF( BUILD_FOR_WEB_DEPLOYMENT )

  IF( BUILD_AS_OSX_APP OR BUILD_AS_ELECTRON_APP OR BUILD_AS_LOCAL_SERVER OR BUILD_AS_UNIT_TEST_SUITE )
    set( sources ${sources} src/SpecFileQuery.cpp src/SpecFileQueryWidget.cpp src/SpecFileQueryDbCache.cpp )
    set( headers ${headers} InterSpec/SpecFileQuery.h InterSpec/SpecFileQueryWidget.h InterSpec/SpecFileQueryDbCache.h )
    list( APPEND sources js/SpecFileQueryWidget.js )

    add_subdirectory( external_libs/pugixml-1.9 )
    list( APPEND LIBRARIES_TO_LINK_TO pugixml )
  ENDIF( BUILD_AS_OSX_APP OR BUILD_AS_ELECTRON_APP OR BUILD_AS_LOCAL_SERVER OR BUILD_AS_UNIT_TEST_SUITE )
ENDIF( USE_SPECRUM_FILE_QUERY_WIDGET )

if( USE_SPECTRUM_CHART_D3 )
   list( APPEND sources
      src/D3SpectrumDisplayDiv.cpp
      external_libs/SpecUtils/d3_resources/d3.v3.min.js
      external_libs/SpecUtils/d3_resources/c.min.js
      external_libs/SpecUtils/d3_resources/SpectrumChartD3.js
      external_libs/SpecUtils/d3_resources/SpectrumChartD3.css
      )
    set( headers ${headers} InterSpec/D3SpectrumDisplayDiv.h )
endif( USE_SPECTRUM_CHART_D3 )


FIND_PATH( WT_RESOURCES_DIRECTORY form.css
           ${Wt_INCLUDE_DIR}/../share/Wt/resources )

set( SpecUtils_REBIN_FILES_TO_SINGLE_BINNING ON CACHE BOOL "Set Rebin spectrum to single energy cal" )
set( SpecUtils_NO_BOOST_LIB OFF CACHE BOOL "Use boost for threading primitives in SpecUtils" )
set( SpecUtils_USE_WT_THREADPOOL ON CACHE BOOL "" )
set( SpecUtils_EXTERNALLY_DEFINED_LOG_MESSAGE ON CACHE BOOL "" )
set( SpecUtils_BUILD_TEST OFF CACHE BOOL "" )
set( SpecUtils_PYTHON_BINDINGS OFF CACHE BOOL "" )
set( SpecUtils_JAVA_SWIG OFF CACHE BOOL "" )
set( SpecUtils_D3_SUPPORT_FILE_STATIC OFF CACHE BOOL "" )
include_directories( "${Wt_INCLUDE_DIR}" )

add_subdirectory( external_libs/SpecUtils )
list( APPEND LIBRARIES_TO_LINK_TO SpecUtils )

add_subdirectory( external_libs/SandiaDecay )
list( APPEND LIBRARIES_TO_LINK_TO SandiaDecay )

set( RAPID_XML_INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/3rdparty )
#find_path( RAPID_XML_INC_DIR rapidxml/rapidxml.hpp ${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/3rdparty )
INCLUDE_DIRECTORIES( ${RAPID_XML_INC_DIR} )

#We should have figured out all the header and source files we need by here, so
#  lets add the executable to the build
if( ANDROID OR IOS )
  if( IOS )
    set( sources ${sources} src/InterSpecServer.cpp )
    set( headers ${headers} InterSpec/InterSpecServer.h )
    add_library( ${LIBRARYNAME} STATIC ${sources} ${headers} )
  else( IOS )
    add_library( ${LIBRARYNAME} SHARED main.cpp ${sources} ${headers} )
    set( EXECUTABLE_NAME ${LIBRARYNAME} )
  endif(IOS)
else( ANDROID OR IOS )
  if( BUILD_AS_OSX_APP )

    set( sources ${sources} src/InterSpecServer.cpp )
    set( headers ${headers} InterSpec/InterSpecServer.h )

    SET(MACOSX_BUNDLE_INFO_STRING "InterSpec - Sandia National Lab, Will Johnson" )
    SET(MACOSX_BUNDLE_ICON_FILE InterSpec )
    SET(MACOSX_BUNDLE_LONG_VERSION_STRING "" )
    SET(MACOSX_BUNDLE_BUNDLE_NAME "InterSpec" )
    SET(MACOSX_BUNDLE_SHORT_VERSION_STRING "InterSpec" )
    SET(MACOSX_BUNDLE_BUNDLE_VERSION "0" )
    SET(MACOSX_BUNDLE_COPYRIGHT "Sandia National Labs, Will Johnson" )
    SET(MACOSX_BUNDLE_GUI_IDENTIFIER "gov.sandia.InterSpec" )
    SET(MACOSX_BUNDLE_INFO_PLIST  )

    SET( ${PRODUCT_NAME} "InterSpec" )

    #${MACOSX_DEPLOYMENT_TARGET}

    SET( MACOSX_BUNDLE_BUNDLE_NAME "InterSpec" )

    FIND_LIBRARY(WEBKIT_LIBRARY WebKit)
    FIND_LIBRARY(COCOA_LIBRARY Cocoa)
    FIND_LIBRARY(COREDATA_LIBRARY CoreData)
    FIND_LIBRARY(APPKIT_LIBRARY AppKit)
    FIND_LIBRARY(FOUNDATION_LIBRARY Foundation)

    MARK_AS_ADVANCED( WEBKIT_LIBRARY COCOA_LIBRARY COREDATA_LIBRARY APPKIT_LIBRARY FOUNDATION_LIBRARY )
    list( APPEND LIBRARIES_TO_LINK_TO
          ${WEBKIT_LIBRARY} ${COCOA_LIBRARY} ${APPKIT_LIBRARY}
          ${COREDATA_LIBRARY} ${FOUNDATION_LIBRARY} )

    set( app_sources
          ${CMAKE_CURRENT_SOURCE_DIR}/target/osx/InterSpec.icns
          "target/osx/en.lproj/MainMenu.xib"
          "target/osx/AppDelegate.mm"
          "target/osx/main.m" )

    set( app_headers ${headers} "target/osx/AppDelegate.h" )

    ADD_EXECUTABLE( ${EXECUTABLE_NAME} MACOSX_BUNDLE ${app_sources} ${app_headers} ${sources} ${headers} )
    SET_TARGET_PROPERTIES( ${EXECUTABLE_NAME} PROPERTIES MACOSX_BUNDLE TRUE )

    SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/target/osx/InterSpec.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

    SET( EN_LOCALIZATION_FILES InfoPlist.strings Credits.rtf )

    #set( LINK_TO_OBJECT "${EXECUTABLE_NAME}" )
    #target_link_libraries( ${EXECUTABLE_NAME} PRIVATE ${LIBRARIES_TO_LINK_TO} )
    #set_property(TARGET ${EXECUTABLE_NAME} PROPERTY CXX_STANDARD 11)
    #set_property(TARGET prog PROPERTY CXX_STANDARD_REQUIRED ON)

    add_custom_target(copy)
    foreach(_file ${EN_LOCALIZATION_FILES})
        add_custom_command(TARGET copy PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
            copy "${CMAKE_CURRENT_SOURCE_DIR}/target/osx/en.lproj/${_file}"
                  "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/\${CONFIGURATION}/${EXECUTABLE_NAME}.app/Contents/Resources/en.lproj/${_file}" )
    endforeach()

    SET( SUPPORT_DIRECTORIES InterSpec_resources data example_spectra external_libs/SpecUtils/d3_resources  )
    foreach(_dir ${SUPPORT_DIRECTORIES})
      add_custom_command(TARGET copy PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_directory ${PROJECT_SOURCE_DIR}/${_dir}
        ${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/\${CONFIGURATION}/${EXECUTABLE_NAME}.app/Contents/Resources/${_dir} )
    endforeach()

    add_custom_command(TARGET copy PRE_BUILD COMMAND ${CMAKE_COMMAND} -E
        copy_directory ${WT_RESOURCES_DIRECTORY}
        ${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/\${CONFIGURATION}/${EXECUTABLE_NAME}.app/Contents/Resources/resources )


    find_program(IBTOOL ibtool HINTS "/usr/bin" "${OSX_DEVELOPER_ROOT}/usr/bin")
    if( ${IBTOOL} STREQUAL "IBTOOL-NOTFOUND" )
        message(SEND_ERROR "could fine ibtool to compile the .xib files")
    endif()

    add_custom_target( makenib )
    add_custom_command( TARGET makenib POST_BUILD
        COMMAND ${IBTOOL} --errors --warnings --notices --output-format human-readable-text
                --compile "${CMAKE_CURRENT_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/\${CONFIGURATION}/${EXECUTABLE_NAME}.app/Contents/Resources/en.lproj/MainMenu.nib"
                    "${CMAKE_CURRENT_SOURCE_DIR}/target/osx/en.lproj/MainMenu.xib"
                COMMENT "Compiling MainMenu.xib" )

    add_dependencies( makenib copy )
    add_dependencies( ${EXECUTABLE_NAME} makenib )

    set_target_properties( ${EXECUTABLE_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST
          "${CMAKE_CURRENT_SOURCE_DIR}/target/osx/Info.plist.template" )

  else( BUILD_AS_OSX_APP )

    if( BUILD_AS_UNIT_TEST_SUITE )
      add_library( ${LIBRARYNAME} STATIC ${sources} ${headers} )
      #list( APPEND LIBRARIES_TO_LINK_TO ${LIBRARYNAME} )
      add_executable( ${EXECUTABLE_NAME} ${GUI_TYPE} main.cpp )
    else( BUILD_AS_UNIT_TEST_SUITE )
      add_executable( ${EXECUTABLE_NAME} ${GUI_TYPE} main.cpp ${sources} ${headers} )
    endif( BUILD_AS_UNIT_TEST_SUITE )

  endif( BUILD_AS_OSX_APP )
endif( ANDROID OR IOS )



if(ANDROID)
    INCLUDE_DIRECTORIES( target )
    set( headers ${headers} target/android/AndroidUtils.hpp )
ENDIF(ANDROID)


add_subdirectory( external_libs/Minuit2 )
list( APPEND LIBRARIES_TO_LINK_TO minuit2 )


if( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE OR BUILD_AS_UNIT_TEST_SUITE OR BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)
    list( APPEND LIBRARIES_TO_LINK_TO ${Wt_TEST_LIBRARY} )
endif( BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE OR BUILD_AS_UNIT_TEST_SUITE OR BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT )


IF( IOS )
  #kludge for Xcode linker error
  set( HAS_WTDBOMYSQL OFF )
ENDIF( IOS )

IF( USE_MYSQL_DB AND USE_SQLITE3_DB )
  MESSAGE(FATAL_ERROR "Only one type of database (MySQL or SQLITE3) may be chosen" )
ENDIF( USE_MYSQL_DB AND USE_SQLITE3_DB )

IF( NOT USE_MYSQL_DB AND NOT USE_SQLITE3_DB )
  MESSAGE(FATAL_ERROR "You must choose either a MySQL or SQLITE3 database" )
ENDIF( NOT USE_MYSQL_DB AND NOT USE_SQLITE3_DB )


INCLUDE_DIRECTORIES(
    ${Wt_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}
    external_libs
    .
)


if( NOT IOS )

IF(WIN32)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    list( APPEND LIBRARIES_TO_LINK_TO ${Wt_DEBUG_LIBRARIES} )
  else(CMAKE_BUILD_TYPE STREQUAL "Debug")
    list( APPEND LIBRARIES_TO_LINK_TO ${Wt_LIBRARIES} )
  endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  link_directories( ${Boost_LIBRARY_DIR} )

  IF( Boost_CHRONO_FOUND )
    list( APPEND LIBRARIES_TO_LINK_TO ${Boost_CHRONO_LIBRARY} )
  ENDIF( Boost_CHRONO_FOUND )

ELSE(WIN32)
  list( APPEND LIBRARIES_TO_LINK_TO ${Wt_LIBRARIES} )
ENDIF(WIN32)

  list( APPEND LIBRARIES_TO_LINK_TO
      ${Boost_LIBRARIES}
      ${Boost_SYSTEM_LIBRARY}
      ${Boost_REGEX_LIBRARY}
      ${Boost_THREAD_LIBRARY}
      ${Boost_FILESYSTEM_LIBRARY}
      ${Boost_RANDOM_LIBRARY}
      ${Boost_PROGRAM_OPTIONS_LIBRARY}
      ${Boost_DATE_TIME_LIBRARY}
  )

  find_package( Threads REQUIRED )
  target_link_libraries( ${LINK_TO_OBJECT} PUBLIC Threads::Threads )

endif( NOT IOS )
#ENDIF(WIN32)

if(APPLE AND NOT IOS)
  set( CMAKE_OSX_DEPLOYMENT_TARGET "10.10" )
endif()

#20180104: on my mac I tried compiling Wt with support for PDF (libharu), but only
#  had the static libhpdfs, so I had to add bellow to link to this
FIND_LIBRARY( HPDF_LIBRARY NAME hpdfs HINTS ${Wt_INCLUDE_DIR}/../lib )
if( NOT HPDF_LIBRARY )
  message( "Could not find libharu - if Wt was built with PDF support, you will get link errors" )
else()
  list( APPEND LIBRARIES_TO_LINK_TO ${HPDF_LIBRARY} )
  FIND_LIBRARY( PNG_LIBRARY NAME png HINTS ${Wt_INCLUDE_DIR}/../lib )
  if( NOT PNG_LIBRARY )
    message( FATAL "Could not find libpng - requred for libharu" )
  else()
    list( APPEND LIBRARIES_TO_LINK_TO ${PNG_LIBRARY} )
  endif()
endif()


#bellow is what I need to get things working 20131120 for statically linked
if( TRY_TO_STATIC_LINK )
  if(APPLE AND NOT IOS)
    set( CMAKE_OSX_DEPLOYMENT_TARGET "10.10" )

    FIND_LIBRARY(LIBZ_LIBRARY z)

    list( APPEND LIBRARIES_TO_LINK_TO
         ${Boost_RANDOM_LIBRARY}
         ${Boost_PROGRAM_OPTIONS_LIBRARY}
         ${LIBZ_LIBRARY}
         /usr/lib/libSystem.B.dylib
         /System/Library/Frameworks/Carbon.framework/Versions/A/Carbon
         /System/Library/Frameworks/Foundation.framework/Versions/C/Foundation
    )

    if(USE_OSX_NATIVE_MENU)
      list( APPEND LIBRARIES_TO_LINK_TO /System/Library/Frameworks/AppKit.framework/Versions/C/AppKit )
    endif(USE_OSX_NATIVE_MENU)
  endif(APPLE AND NOT IOS)

  list( APPEND LIBRARIES_TO_LINK_TO ${CMAKE_DL_LIBS} )

  if( USE_MYSQL_DB )
    include( cmake/WtFindMysql.txt )
    list( APPEND LIBRARIES_TO_LINK_TO ${MYSQL_LIBRARIES} )
  endif( USE_MYSQL_DB )

  if( Wt_FCGI_LIBRARY )
    include( cmake/WtFindFcgi.txt )
    list( APPEND LIBRARIES_TO_LINK_TO ${FCGI_LIBRARIES} )
  endif( Wt_FCGI_LIBRARY )

  include( cmake/WtFindSsl.txt )
  if( SSL_FOUND )
    list( APPEND LIBRARIES_TO_LINK_TO ${SSL_LIBRARIES} )
  endif( SSL_FOUND )

  #link to zlib since Wt might
  if("${CMAKE_SYSTEM}" MATCHES "Linux" AND NOT SUPPORT_ZIPPED_SPECTRUM_FILES)
    list( APPEND LIBRARIES_TO_LINK_TO ${ZLIB_LIBRARIES} )
  endif("${CMAKE_SYSTEM}" MATCHES "Linux" AND NOT SUPPORT_ZIPPED_SPECTRUM_FILES)

  if("${CMAKE_SYSTEM}" MATCHES "Linux" AND NOT ANDROID)
    list( APPEND LIBRARIES_TO_LINK_TO -lrt )
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++" )
  endif("${CMAKE_SYSTEM}" MATCHES "Linux" AND NOT ANDROID)
endif( TRY_TO_STATIC_LINK )


IF( BUILD_AS_UNIT_TEST_SUITE )
  enable_testing()
  list( APPEND LIBRARIES_TO_LINK_TO ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
ENDIF( BUILD_AS_UNIT_TEST_SUITE )


IF( PERFORM_DEVELOPER_CHECKS )
    list( APPEND LIBRARIES_TO_LINK_TO ${Boost_DATE_TIME_LIBRARY} )
ENDIF( PERFORM_DEVELOPER_CHECKS )


set(CMAKE_REQUIRED_INCLUDES ${Wt_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
INCLUDE(CheckIncludeFileCXX)

if( Wt_DBOSQLITE3_LIBRARY )
  check_include_file_cxx( "Wt/Dbo/backend/Sqlite3"  HAS_WTDBOSQLITE3 )
endif( Wt_DBOSQLITE3_LIBRARY )

if( Wt_MYSQL_LIBRARY )
  check_include_file_cxx( "Wt/Dbo/backend/MySQL"    HAS_WTDBOMYSQL )
endif( Wt_MYSQL_LIBRARY )

if( Wt_DBOPOSTGRES_LIBRARY )
  check_include_file_cxx( "Wt/Dbo/backend/Postgres" HAS_WTDBOPOSTGRES )
endif( Wt_DBOPOSTGRES_LIBRARY )

if( Wt_FIREBIRD_LIBRARY )
  check_include_file_cxx( "Wt/Dbo/backend/Firebird" HAS_WTDBOFIREBIRD )
endif( Wt_FIREBIRD_LIBRARY )

INCLUDE( CheckIncludeFiles )

if( ZLIB_FOUND AND Boost_IOSTREAMS_FOUND )
  #CHECK_INCLUDE_FILES( "zlib.h" HAS_ZLIB_SUPPORT )
  list( APPEND LIBRARIES_TO_LINK_TO ${ZLIB_LIBRARIES} ${Boost_IOSTREAMS_LIBRARY} )
endif( ZLIB_FOUND AND Boost_IOSTREAMS_FOUND )

if( SUPPORT_ZIPPED_SPECTRUM_FILES OR ENABLE_RESOURCE_UPDATES )
#CMake seems to fail to find zlib on ios/android (even using find_package(ZLIB REQUIRED)), but we know it will be in
#  the systems path anyway, so just blindly link to zlib
  if( NOT (IOS OR ANDROID) )
    if( NOT ZLIB_FOUND )
      MESSAGE(FATAL_ERROR "ZLIB required to enable SUPPORT_ZIPPED_SPECTRUM_FILES or ENABLE_RESOURCE_UPDATES, but it wasnt found" )
    else( NOT ZLIB_FOUND )
      include_directories( ${ZLIB_INCLUDE_DIRS} )
      list( APPEND LIBRARIES_TO_LINK_TO ${ZLIB_LIBRARIES} )
    endif( NOT ZLIB_FOUND )
  else( NOT (IOS OR ANDROID) )
    list( APPEND LIBRARIES_TO_LINK_TO z )
  endif( NOT (IOS OR ANDROID) )
endif( SUPPORT_ZIPPED_SPECTRUM_FILES OR ENABLE_RESOURCE_UPDATES )


IF( Wt_MYSQL_LIBRARY AND BUILD_FOR_WEB_DEPLOYMENT AND USE_MYSQL_DB )
  set( DATABASE_PASSWORD_FILE "/path/to/passwords/databases.xml" CACHE STRING
         "Location of the XML passwords file for connecting to MySQL database" )
ENDIF( Wt_MYSQL_LIBRARY AND BUILD_FOR_WEB_DEPLOYMENT AND USE_MYSQL_DB )

IF( USE_MYSQL_DB AND NOT HAS_WTDBOMYSQL )
  MESSAGE(FATAL_ERROR "You do not have the MySQL Wt Dbo backend installed, but you have asked to use it")
ENDIF( USE_MYSQL_DB AND NOT HAS_WTDBOMYSQL )

SET( HAS_WTDBOSQLITE3 1 )
IF( USE_SQLITE3_DB AND NOT HAS_WTDBOSQLITE3 )
  MESSAGE(FATAL_ERROR "You do not have the SQLITE3 Wt Dbo backend installed, but you have asked to use it")
ENDIF( USE_SQLITE3_DB AND NOT HAS_WTDBOSQLITE3 )


configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec/InterSpec_config.h.in
    ${CMAKE_BINARY_DIR}/InterSpec_config.h
)

add_subdirectory( external_libs/Cuba-3.0 )
list( APPEND LIBRARIES_TO_LINK_TO Cuba-3.0 )

add_subdirectory( external_libs/muparserx-4.0.7 )
include_directories( external_libs/muparserx-4.0.7/parser )
list( APPEND LIBRARIES_TO_LINK_TO muparserx )

IF(ANDROID)
  list( APPEND LIBRARIES_TO_LINK_TO log )
ENDIF(ANDROID)


message( "Linking ${LINK_TO_OBJECT} to ${LIBRARIES_TO_LINK_TO}" )
target_link_libraries( ${LINK_TO_OBJECT} PUBLIC ${LIBRARIES_TO_LINK_TO} )


IF( BUILD_AS_UNIT_TEST_SUITE )
  target_link_libraries( ${EXECUTABLE_NAME} PUBLIC ${LIBRARYNAME} )

  add_executable( testCalibration.exe testing/testCalibration.cpp )
  target_link_libraries( testCalibration.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
  add_test( "\"Test Calibration\""
                 ${EXECUTABLE_OUTPUT_PATH}/testCalibration.exe
                 "--indir=${testdir}" --log_level=test_suite --run_test=testCalibration --catch_system_error=yes ) 
  add_test( "\"Test Calibration\""
                 ${EXECUTABLE_OUTPUT_PATH}/testCalibration.exe
                 "--indir=${testdir}" --log_level=test_suite --run_test=testFindEnergy --catch_system_error=yes )

  add_executable( testFileOpen.exe testing/testFileOpen.cpp )
  target_link_libraries( testFileOpen.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )


  SET( curdir ${PROJECT_SOURCE_DIR}/testing/interspec_test_spectra )
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  FOREACH(child ${children})
    SET( testdir ${curdir}/${child} )
    IF(IS_DIRECTORY ${testdir})
       #uhg, if we dont get rid of spaces in the file name, it gets interprested as a
       #  seperate command line argument - even if we escape it with quotes
       string(REPLACE " " "%20" testdir ${testdir} )
       message( ${testdir} )
       add_test( "\"Test open ${child}\""
                 ${EXECUTABLE_OUTPUT_PATH}/testFileOpen.exe 
                 "--indir=${testdir}" --log_level=test_suite --run_test=testFileOpen --catch_system_error=yes ) 
    ENDIF()
  ENDFOREACH()


add_executable( testIsoTimeString.exe testing/testIsoTimeString.cpp )
  target_link_libraries( testIsoTimeString.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
add_test( TestIsoTimeString ${EXECUTABLE_OUTPUT_PATH}/testIsoTimeString.exe "--indir=${testdir}" --log_level=test_suite --run_test=testIsoTimeString --catch_system_error=yes )

add_executable( testSafeGetLine.exe testing/testSafeGetLine.cpp )
  target_link_libraries( testSafeGetLine.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
add_test( TestSafeGetLine ${EXECUTABLE_OUTPUT_PATH}/testSafeGetLine.exe "--indir=${testdir}" --log_level=test_suite --run_test=testSafeGetLine --catch_system_error=yes )

add_executable( testTimeFromString.exe testing/testTimeFromString.cpp )
  target_link_libraries( testTimeFromString.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
add_test( TestTimeFromString ${EXECUTABLE_OUTPUT_PATH}/testTimeFromString.exe "--indir=${testdir}" --log_level=test_suite --run_test=testTimeFromString --catch_system_error=yes )


add_executable( testUtilityStringFunctions.exe testing/testUtilityStringFunctions.cpp )
  target_link_libraries( testUtilityStringFunctions.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
add_test( "\"Test Utility String Functions\"" ${EXECUTABLE_OUTPUT_PATH}/testUtilityStringFunctions.exe "--indir=${testdir}" --log_level=test_suite --run_test=testUtilityStringFunctions --catch_system_error=yes )

add_test( "\"Test Utility String Functions\"" ${EXECUTABLE_OUTPUT_PATH}/testUtilityStringFunctions.exe "--indir=${testdir}" --log_level=test_suite --run_test=testUtilityStringFunctions --catch_system_error=yes ) 


add_executable( testPhysicalUnits.exe testing/testPhysicalUnits.cpp )
target_link_libraries( testPhysicalUnits.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
add_test( "\"Test PhysicalUnits Functions\"" ${EXECUTABLE_OUTPUT_PATH}/testPhysicalUnits.exe "--indir=${testdir}" --log_level=test_suite --run_test=testPhysicalUnits_suite --catch_system_error=yes )

add_test( "\"Test Utility String Functions\"" ${EXECUTABLE_OUTPUT_PATH}/testUtilityStringFunctions.exe "--indir=${testdir}" --log_level=test_suite --run_test=testUtilityStringFunctions --catch_system_error=yes )

  
add_executable( test_split_to_floats_and_ints.exe testing/test_split_to_floats_and_ints.cpp )
  target_link_libraries( test_split_to_floats_and_ints.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )


add_executable( test_utf8_str_size_limit.exe testing/test_utf8_str_size_limit.cpp )
  target_link_libraries( test_utf8_str_size_limit.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )

add_executable( test_utf8_limit_str_size.exe testing/test_utf8_limit_str_size.cpp )
  target_link_libraries( test_utf8_limit_str_size.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )

add_executable( test_utf8_str_len.exe testing/test_utf8_str_len.cpp )
  target_link_libraries( test_utf8_str_len.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
 

 add_test( TestSplittoFloats ${EXECUTABLE_OUTPUT_PATH}/test_split_to_floats_and_ints.exe "--indir=${testdir}" --log_level=test_suite --run_test=split_to_floats_and_ints_suite --catch_system_error=yes )

add_test( TestUtf8LimitStrSize ${EXECUTABLE_OUTPUT_PATH}/test_utf8_limit_str_size.exe "--indir=${testdir}" --log_level=test_suite --run_test=test_utf8_limit_str_size --catch_system_error=yes ) 

add_test( TestUtf8StrSizeLimit ${EXECUTABLE_OUTPUT_PATH}/test_utf8_str_size_limit.exe "--indir=${testdir}" --log_level=test_suite --run_test=test_utf8_str_size_limit --catch_system_error=yes ) 

add_test( TestUtf8StrLen ${EXECUTABLE_OUTPUT_PATH}/test_utf8_str_len.exe "--indir=${testdir}" --log_level=test_suite --run_test=test_utf8_str_len --catch_system_error=yes )

  add_executable( testDetectorResponseOpen.exe testing/testDetectorResponseOpen.cpp )
  target_link_libraries( testDetectorResponseOpen.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )

  SET( curdir ${PROJECT_SOURCE_DIR}/data/OUO_GadrasDetectors )
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*/Detector.dat)
  FOREACH(child ${children})
    SET( curfile "${curdir}/${child}" )
    string(REPLACE " " "%20" curfile ${curfile} )
    add_test( "\"Detector Response Open ${child}\""
              "${EXECUTABLE_OUTPUT_PATH}/testDetectorResponseOpen.exe"
              ARGS --infile=${curfile} --log_level=test_suite --run_test=testDetectorResponseOpen --catch_system_error=yes )
  ENDFOREACH()
  
  
  add_executable( testDetectorResponseDoesntOpen.exe testing/testDetectorResponseDoesntOpen.cpp )
  target_link_libraries( testDetectorResponseDoesntOpen.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )

  SET( curdir ${PROJECT_SOURCE_DIR}/testing/detector_response_files/shouldnt_be_valid )
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  FOREACH(child ${children})
    SET( curfile "${curdir}/${child}" )
    string(REPLACE " " "%20" curfile ${curfile} )
    message( "Got child ${curfile}" )
    add_test( "\"Detector Response Shouldnt Open ${child}\""
              "${EXECUTABLE_OUTPUT_PATH}/testDetectorResponseDoesntOpen.exe"
              ARGS --infile=${curfile} --log_level=test_suite --run_test=testDetectorResponseDoesntOpen --catch_system_error=yes )
  ENDFOREACH()

  add_executable( testUtilityFilesystemFunctions.exe testing/testUtilityFilesystemFunctions.cpp )
  target_link_libraries( testUtilityFilesystemFunctions.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
  add_test( testUtilityFilesystemFunctions ${EXECUTABLE_OUTPUT_PATH}/testUtilityFilesystemFunctions.exe  --log_level=test_suite --run_test=testUtilityFilesystemFunctions --catch_system_error=yes )

#  add_executable( peakFitCompare.exe testing/peakFitCompare.cpp )
#  target_link_libraries( peakFitCompare.exe PRIVATE ${LIBRARIES_TO_LINK_TO} ${LIBRARYNAME} )
#  add_test( "\"Test peak fitting\""
#              "${EXECUTABLE_OUTPUT_PATH}/peakFitCompare.exe"
#              ARGS --log_level=test_suite --run_test=peakFitCompare --catch_system_error=yes )


ENDIF( BUILD_AS_UNIT_TEST_SUITE )


set( OTHER_SUPPORT_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/README.md
     ${CMAKE_CURRENT_SOURCE_DIR}/data/default_preferences.xml
     ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/themes/dark/dark.css
     ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/copyright_and_about.xml
     ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/use_instructions.xml
  )

  FILE(GLOB CSS_FILES ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/*.css )
  list( APPEND OTHER_SUPPORT_FILES ${CSS_FILES} )

  FILE(GLOB HELP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/*help.xml )
  list( APPEND OTHER_SUPPORT_FILES ${HELP_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/help.json )

if( USE_SIMPLE_NUCLIDE_ASSIST )
  list( APPEND OTHER_SUPPORT_FILES "InterSpec_resources/SimpleNuclideAssist.css" )
endif( USE_SIMPLE_NUCLIDE_ASSIST )


if( USE_SEARCH_MODE_3D_CHART )
  list( APPEND OTHER_SUPPORT_FILES "InterSpec_resources/SearchMode3DChart.css" )
endif( USE_SEARCH_MODE_3D_CHART )


add_custom_target( interspec_resources SOURCES ${OTHER_SUPPORT_FILES} )


IF( BUILD_AS_ELECTRON_APP )

  find_path( WEBSOCKETPP_INCLUDE_DIR NAMES websocketpp/version.hpp REQUIRED
             HINTS "${CMAKE_CURRENT_SOURCE_DIR}/external_libs/websocketpp/" )

  include_directories( ${WEBSOCKETPP_INCLUDE_DIR} )

  set( ELECTRON_SUPPORT_FILES
       ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/main.js
       ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/loading.html
       ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/README.md
       ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/package.json
       ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/copy_resources.in.js
       ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/launch_options.json
      )

  add_custom_target( electron_resources SOURCES ${ELECTRON_SUPPORT_FILES} )

#We cant get the executable location at configure time, so I'll just hardcode
#  this for now, for a workaround see
#  https://cmake.org/pipermail/cmake/2012-May/050221.html
#  (It would be nice to use $<TARGET_FILE:${EXECUTABLE_NAME}>)
#See the script : target/electron/app/copy_resources.in.js for use of this variabel
  IF(WIN32)
      set( INTERSPEC_EXE_LOCATION "${PROJECT_BINARY_DIR}/MinSizeRel/InterSpec.exe.exe")
  ELSE(WIN32)
       set( INTERSPEC_EXE_LOCATION "${PROJECT_BINARY_DIR}/bin/InterSpec.exe")
  ENDIF(WIN32)

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/copy_resources.in.js
    ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/copy_resources.js
  )

ENDIF( BUILD_AS_ELECTRON_APP )

#Copy files to the install directory
set( OTHER_FILES runLocalOn8080.sh )

foreach(_file ${OTHER_FILES})
  configure_file( ${PROJECT_SOURCE_DIR}/${_file} ${PROJECT_BINARY_DIR}/${_file} COPYONLY )
endforeach()


IF( IOS )
  MESSAGE("CMAKE_CURRENT_SOURCE_DIR is ${CMAKE_CURRENT_SOURCE_DIR}")

  MESSAGE("Looking for form.css in ${CMAKE_CURRENT_SOURCE_DIR}/wt-3.3.1/resources")
  IF( NOT WT_RESOURCES_DIRECTORY )
      FIND_PATH( WT_RESOURCES_DIRECTORY form.css HINTS ${CMAKE_CURRENT_SOURCE_DIR}/wt-3.3.1/resources NO_DEFAULT_PATHS NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_FIND_ROOT_PATH )
  ENDIF( NOT WT_RESOURCES_DIRECTORY )
  MESSAGE("***** WT_RESOURCES_DIRECTORY is ${WT_RESOURCES_DIRECTORY}")

  MESSAGE("WT_RESOURCES_DIRECTORY IS ${WT_RESOURCES_DIRECTORY}")

#next line doesnt seem to work?
#We would also want to do this for Wt and and the actual InterSpec.app projects.
  set_xcode_property( ${LINK_TO_OBJECT} IPHONEOS_DEPLOYMENT_TARGET "7.0" )
ENDIF( IOS )


IF(NOT WT_RESOURCES_DIRECTORY)
  MESSAGE("Could not find the Wt resources path to copy into current directory")
ELSE(NOT WT_RESOURCES_DIRECTORY)
  MESSAGE("Using Wt resources from ${WT_RESOURCES_DIRECTORY}")
  file(COPY ${WT_RESOURCES_DIRECTORY}
      DESTINATION ${PROJECT_BINARY_DIR}/
      )
ENDIF(NOT WT_RESOURCES_DIRECTORY)

IF( NOT ANDROID )
  SET( SUPPORT_DIRECTORIES InterSpec_resources data example_spectra external_libs/SpecUtils/d3_resources )
endif( NOT ANDROID )

foreach(_dir ${SUPPORT_DIRECTORIES})
    IF(WIN32)
#        file(COPY ${PROJECT_SOURCE_DIR}/${_dir} DESTINATION ${PROJECT_BINARY_DIR}/ )
    ELSE(WIN32)
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        ${PROJECT_SOURCE_DIR}/${_dir}
                        ${PROJECT_BINARY_DIR}/${_dir} )
    ENDIF(WIN32)
endforeach()

IF( INCLUDE_ANALYSIS_TEST_SUITE )
    IF(WIN32)
#        file(COPY ${PROJECT_SOURCE_DIR}/testing/analysis_tests DESTINATION ${PROJECT_BINARY_DIR}/analysis_tests )
    ELSE(WIN32)
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        ${PROJECT_SOURCE_DIR}/testing/analysis_tests
                        ${PROJECT_BINARY_DIR}/analysis_tests )
    ENDIF(WIN32)
ENDIF()

configure_file( ${PROJECT_SOURCE_DIR}//external_libs/SandiaDecay/sandia.decay.nocoinc.min.xml ${PROJECT_BINARY_DIR}/data/sandia.decay.xml COPYONLY )

IF( BUILD_AS_UNIT_TEST_SUITE )
    IF(WIN32)
#        file(COPY ${PROJECT_SOURCE_DIR}/testing/interspec_test_spectra
 #                 ${PROJECT_BINARY_DIR}/interspec_test_spectra )
    ELSE(WIN32)
        execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        ${PROJECT_SOURCE_DIR}/testing/interspec_test_spectra
                        ${PROJECT_BINARY_DIR}/interspec_test_spectra )
    ENDIF(WIN32)
ENDIF()

