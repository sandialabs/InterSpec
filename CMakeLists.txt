# This CMakeLists potentially defines two items, InterSpecLib and for all targets
# excpet iOS, Android, and Electron, InterSpecExe

cmake_policy(SET CMP0048 NEW)
project(InterSpec VERSION 1.0.8)

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
set(CMAKE_COLOR_MAKEFILE ON CACHE BOOL "Color" FORCE)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32)
    set(MSVC_RUNTIME "static")
    set(TRY_TO_STATIC_LINK "ON")
    include(cmake/ConfigureMsvc.txt)
    configure_msvc_runtime()

    set(BOOST_PREFIX "C:/install/msvc2017/x64/boost_1_65_1")
    set(ZLIB_PREFIX "C:/install/msvc2017/x64/zlib")
    list(APPEND PATH_SUFFIXES "C:/install/msvc2017/x64/wt-3.3.4")

    set(EXECUTABLE_OUTPUT_PATH . CACHE PATH "Path to executables" FORCE)
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} /bigobj -D_SCL_SECURE_NO_WARNINGS /MP /wd4996 /wd4267 /DWINVER=0x0601 /D_WIN32_WINNT=0x0601"
    )
    #0x0601==Win7, 0x0501==WinXP
else(WIN32)
    set(EXECUTABLE_OUTPUT_PATH bin CACHE PATH "Path to executables" FORCE)
endif(WIN32)

SET(MYSQL_DATABASE_TO_USE "dev" CACHE STRING "Which mysql database to use (if applicable): dev/qc/prod")

option(BUILD_AS_OSX_APP "Build a OSX native app" OFF)
option(
    BUILD_AS_ELECTRON_APP
    "Configures the executable to be ran as a Electron app"
    OFF
)
option(
    BUILD_FOR_WEB_DEPLOYMENT
    "Selects options appropriate for the web in the future"
    OFF
)
option(
    BUILD_AS_LOCAL_SERVER
    "Build for use as a server for local computer (e.g. run from the command line and connect to with a browser on LOCALHOST)"
    ON
)
option(
    BUILD_AS_UNIT_TEST_SUITE
    "Builds unit tests as well the analysis tests (aka end-to-end tests)"
    OFF
)
option(
    BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE
    "Compiles so executable only does offline testing of the user test states in the database (e.g. end-to-end testing)"
    OFF
)
set(TEST_SUITE_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/testing"
    CACHE STRING
    "Path to directory that contains the \"analysis_tests\" directory for saving N42 test states.  Leave empty for CWD."
)

option(
    BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT
    "Build executable for use while testing new code, not running InterSpec"
    OFF
)
option(
    INCLUDE_ANALYSIS_TEST_SUITE
    "Allow whether user can save and load test spectra"
    OFF
)
option(
    PERFORM_DEVELOPER_CHECKS
    "Performs additional computationally expensive tests during execution"
    off
)
option(
    USE_OSX_NATIVE_MENU
    "Mirrors WMenus with OSX native NSMenu (not fully implemented)"
    off
)
option(
    USE_ELECTRON_NATIVE_MENU
    "Mirrors WMenus with Electrons native Menu (not fully implemented); if off will use HTML based menus"
    OFF
)
option(
    USE_HIGH_BANDWIDTH_INTERACTIONS
    "Allow more interactiveness with the charts"
    ON
)
option(
    SUPPORT_ZIPPED_SPECTRUM_FILES
    "(Work in progress) Allows opening zipped spectrum files"
    ON
)
option(TRY_TO_STATIC_LINK "Try to link to static libs whenever possible" OFF)
option(
    DRAW_GAMMA_LINES_LOG_AND_LIN
    "Have gamma lines indicators be drawn in log mode when y-axis is log"
    OFF
)
option(
    ALLOW_URL_TO_FILESYSTEM_MAP
    "Allows a file ID to be specified in the URL to open a file from the filesystem"
    ON
)
option(USE_DB_TO_STORE_SPECTRA "Use the database to store spectra" ON)
option(
    RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE
    "Render photopeaks server-side rather than client side"
    ON
)
option(
    USE_SPECRUM_FILE_QUERY_WIDGET
    "Enables compilation of the spectrum file query tool (not for web deployment)"
    OFF
)
option(
    USE_FEATURE_MARKER_WIDGET
    "Use the feature marker widget to display feature marker options"
    ON
)
option(
    USE_SEARCH_MODE_3D_CHART
    "Enable 3D chart or searchmode/passthrough data.  Requires Wt >=3.3.4 compiled with WebGL/OpenGL support"
    ON
)
option(USE_TERMINAL_WIDGET "Enables use of terminal widget" ON)
set(SpecUtils_ENABLE_D3_CHART ON
    CACHE BOOL "Enables use of exporting D3 HTML files"
)
option(
    USE_SPECTRUM_CHART_D3
    "(experimental) Replaces Wt SpectrumChart with D3.js chart"
    ON
)
option(
    USE_SIMPLE_NUCLIDE_ASSIST
    "Enable using of Mike E. gamma rays of interest to help with nuclide ID"
    OFF
)

option(
    USE_SQLITE3_DB
    "Use SQLITE3 database to store user preferences and spectra"
    ON
)
option(
    USE_MYSQL_DB
    "Use MySQL database to store user preferences and spectra"
    OFF
)
option(USE_GOOGLE_MAP "Use google maps widget" ON)
option(
    DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK
    "Include support for downloading images of the displayed page"
    OFF
)
option(
    USE_CSS_FLEX_LAYOUT
    "Use CSS Flex/Grid Layout in-place of Wts layouts were possible (work in progress)"
    OFF
)

set(MAX_SPECTRUM_MEMMORY_SIZE_MB 256
    CACHE STRING
    "Amount of memory to allow spectra to take up before trying to offload them onto disk when not in use"
)

set(GOOGLE_MAPS_KEY "" CACHE STRING "Google maps api key.")

set(GUI_TYPE "")

if(WIN32 OR BUILD_AS_OSX_APP OR BUILD_AS_ELECTRON_APP)
    set(TRY_TO_STATIC_LINK "ON")
endif(WIN32 OR BUILD_AS_OSX_APP OR BUILD_AS_ELECTRON_APP)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    #Default to debug
    SET(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )

    # Set the possible values of build type for cmake-gui
    set_property(
        CACHE CMAKE_BUILD_TYPE
        PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo"
    )
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)

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

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    #elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    #elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    #get rid of a million warnings by Wt/Dbo/ptr_impl.h in VS2015
    add_definitions("/wd4297")
    add_definitions("/wd4244")
    #set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "/ignore:4099")
endif()

if(ANDROID)
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
    set(INCLUDE_ANALYSIS_TEST_SUITE OFF)
    set(SHARED_LIBS OFF)
    add_definitions(-DANDROID -DWT_NO_STD_WSTRING=ON -DWT_NO_STD_LOCALE=ON)

    #This is just painful.  I cant get the find_* macros to work with the CMake
    # included as part of Android Studio... After wasting what seems like an entire
    # day, its time for the brute force solution
    if(NOT BOOST_WT_PREFIX_BASE)
        message(
            fatal
            "You must define BOOST_WT_PREFIX_BASE when building Android build"
        )
    endif(NOT BOOST_WT_PREFIX_BASE)
    set(USRPREFIX "${BOOST_WT_PREFIX_BASE}/${ANDROID_ABI}")
    set(BOOST_LIB_START "${USRPREFIX}/lib/libboost_")
    set(BOOST_LIB_END "-clang-darwin-mt-1_65_1.a")

    list(APPEND CMAKE_FIND_ROOT_PATH "${USRPREFIX}")

    set(BOOST_PREFIX "${USRPREFIX}")
    set(Boost_LIBRARY_DIR "${USRPREFIX}/lib")
    set(Boost_INCLUDE_DIR "${USRPREFIX}/include")
    set(Boost_ATOMIC_LIBRARY_RELEASE ${BOOST_LIB_START}atomic${BOOST_LIB_END})
    set(Boost_CHRONO_LIBRARY_RELEASE ${BOOST_LIB_START}chrono${BOOST_LIB_END})
    set(Boost_DATE_TIME_LIBRARY_RELEASE
        ${BOOST_LIB_START}date_time${BOOST_LIB_END}
    )
    set(Boost_FILESYSTEM_LIBRARY_RELEASE
        ${BOOST_LIB_START}filesystem${BOOST_LIB_END}
    )
    set(Boost_IOSTREAMS_LIBRARY_RELEASE
        ${BOOST_LIB_START}iostreams${BOOST_LIB_END}
    )
    set(Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE
        ${BOOST_LIB_START}program_options${BOOST_LIB_END}
    )
    set(Boost_RANDOM_LIBRARY_RELEASE ${BOOST_LIB_START}random${BOOST_LIB_END})
    set(Boost_REGEX_LIBRARY_RELEASE ${BOOST_LIB_START}regex${BOOST_LIB_END})
    set(Boost_SYSTEM_LIBRARY_RELEASE ${BOOST_LIB_START}system${BOOST_LIB_END})
    set(Boost_THREAD_LIBRARY_RELEASE ${BOOST_LIB_START}thread${BOOST_LIB_END})
    set(Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE
        ${BOOST_LIB_START}unit_test_framework${BOOST_LIB_END}
    )
    set(Boost_SIGNALS_LIBRARY_RELEASE ${BOOST_LIB_START}signals${BOOST_LIB_END})
    set(Wt_INCLUDE_DIR=${USRPREFIX}/include)
    set(Wt_LIBRARY=${USRPREFIX}/lib/libwt.a)
    set(Wt_TEST_LIBRARY=${USRPREFIX}/lib/libwttest.a)
    set(Wt_HTTP_LIBRARY=${USRPREFIX}/lib/libwthttp.a)
    set(Wt_DBO_LIBRARY=${USRPREFIX}/lib/libwtdbo.a)
    set(Wt_DBOSQLITE3_LIBRARY=${USRPREFIX}/lib/libwtdbosqlite3.a)
endif(ANDROID)

if(IOS)
    set(CMAKE_XCODE_ATTRIBUTE_IPHONEOS_DEPLOYMENT_TARGET "10.2"
        CACHE STRING "Set the minimum deployment target value."
        FORCE
    )

    #For iOS we have to use libtool instead of ar so that the platform (armv7, armv7s, arm64, etc) are kept track of in the archive
    SET(CMAKE_AR "libtool" )

    #CMake adds a cq after the archive command by default, we need to change this to "-static -o"
    SET(CMAKE_CXX_ARCHIVE_CREATE   "<CMAKE_AR> -static -o <TARGET> <LINK_FLAGS> <OBJECTS>")
    SET(CMAKE_C_ARCHIVE_CREATE   "<CMAKE_AR> -static -o <TARGET> <LINK_FLAGS> <OBJECTS>")

    set(BUILD_AS_LOCAL_SERVER OFF)
    set(ALLOW_URL_TO_FILESYSTEM_MAP ON)
    set(INCLUDE_ANALYSIS_TEST_SUITE OFF)
    set(PERFORM_DEVELOPER_CHECKS OFF)
    set(TRY_TO_STATIC_LINK ON)
    set(USE_SPECRUM_FILE_QUERY_WIDGET OFF)
    set(MAX_SPECTRUM_MEMMORY_SIZE_MB "32")
    set(FRAMEWORKDIR "${CMAKE_CURRENT_SOURCE_DIR}/target/ios")
    set(CMAKE_SYSTEM_FRAMEWORK_PATH "${FRAMEWORKDIR}")
    set(USE_BOOST_FRAMEWORK OFF)
    set(DEFAULT_WT_BOOST_DISCOVERY TRUE)
    set(CMAKE_FRAMEWORK_PATH "${FRAMEWORKDIR}")
    set(BOOST_FRAMEWORK_PATH "${CMAKE_FRAMEWORK_PATH}/framework/boost.framework")
    set(BOOST_PREFIX "${FRAMEWORKDIR}/prefix")
    set(Boost_DIR "${FRAMEWORKDIR}/prefix")
    set(Boost_INCLUDE_DIR "${Boost_DIR}/include")
    set(Boost_LIBRARY_DIRS "${Boost_DIR}/lib")
    set(BOOST_LIBRARYDIR "${Boost_DIR}/lib")
    set(BOOST_ROOT "${Boost_DIR}")
    SET(CMAKE_MODULE_PATH ${Boost_DIR} ${CMAKE_MODULE_PATH} )
    set(SHARED_LIBS OFF)
    set(CMAKE_FIND_ROOT_PATH ${CMAKE_FIND_ROOT_PATH} ${Boost_DIR})
    set(SUPPORT_ZIPPED_SPECTRUM_FILES ON)
    set(SpecUtils_ENABLE_D3_CHART ON)
    set(USE_SPECTRUM_CHART_D3 OFF)
    set(USE_GOOGLE_MAP OFF)

    #if you change this IPHONEOS_DEPLOYMENT_TARGET also change the variable in
    # iosbuild.rb everywhere, and in the Xcode project.
    set(ENV{IPHONEOS_DEPLOYMENT_TARGET} 10.2)

    #We know where Wt is for iOS, so just force it
    if(NOT Wt_INCLUDE_DIR)
        SET(Wt_DIR "${CMAKE_CURRENT_SOURCE_DIR}/target/ios/prefix")
        SET(Wt_INCLUDE_DIR "${Wt_DIR}/include")
    endif(NOT Wt_INCLUDE_DIR)
endif(IOS)

if(BUILD_AS_OSX_APP)
    set(ALLOW_URL_TO_FILESYSTEM_MAP ON)
    set(USE_OSX_NATIVE_MENU ON)
    set(USE_SPECRUM_FILE_QUERY_WIDGET ON)
    set(SpecUtils_ENABLE_D3_CHART ON)
    set(USE_SPECTRUM_CHART_D3 ON)
endif(BUILD_AS_OSX_APP)

if(BUILD_AS_ELECTRON_APP)
   option(
    INTERSPEC_LIBRARY_STATIC
    "Builds the InterSpec library as a static library that will be statically linked into the electron module; disable this if you want to use the shared library to run InterSpec server from C, or Python, or something - this is a little bit of a hack for the moment, and a seperate target should probably be defined for this use-case."
    ON
    )
    
    set(CMAKE_POSITION_INDEPENDENT_CODE ON
        CACHE BOOL "Set PIC to make electron lib"
    )
    set(BUILD_AS_LOCAL_SERVER OFF)
    set(ALLOW_URL_TO_FILESYSTEM_MAP ON)
    set(USE_OSX_NATIVE_MENU OFF)
    set(USE_SPECRUM_FILE_QUERY_WIDGET ON)
    set(SpecUtils_ENABLE_D3_CHART ON)
    set(USE_SPECTRUM_CHART_D3 ON)
    set(INCLUDE_ANALYSIS_TEST_SUITE OFF)
endif(BUILD_AS_ELECTRON_APP)

if(BUILD_AS_UNIT_TEST_SUITE)
    set(BUILD_AS_LOCAL_SERVER OFF)
    set(USE_ENERGY_VS_TIME_PLOTS OFF)
    set(USE_DB_TO_STORE_SPECTRA ON)
    set(ALLOW_URL_TO_FILESYSTEM_MAP ON)
    set(INCLUDE_ANALYSIS_TEST_SUITE ON)
    set(TRY_TO_STATIC_LINK ON)
    set(PERFORM_DEVELOPER_CHECKS ON)
    set(USE_SPECRUM_FILE_QUERY_WIDGET ON)
    set(SpecUtils_ENABLE_D3_CHART ON)
    set(USE_SPECTRUM_CHART_D3 ON)
endif(BUILD_AS_UNIT_TEST_SUITE)

if(BUILD_FOR_WEB_DEPLOYMENT)
    set(TRY_TO_STATIC_LINK ON)
    set(CMAKE_BUILD_TYPE:STRING "Release")
    set(ALLOW_URL_TO_FILESYSTEM_MAP OFF)
    set(USE_DB_TO_STORE_SPECTRA ON)
    #    set( USE_SQLITE3_DB ON )
    #    set( USE_MYSQL_DB OFF )
    set(DATABASE_PASSWORD_FILE "/path/to/passwords/databases.xml")
    set(USE_GOOGLE_MAP OFF)
    set(INCLUDE_ANALYSIS_TEST_SUITE OFF)
    set(PERFORM_DEVELOPER_CHECKS OFF)
    set(SUPPORT_ZIPPED_SPECTRUM_FILES ON)
    set(USE_SPECRUM_FILE_QUERY_WIDGET OFF)
    set(USE_TERMINAL_WIDGET OFF)
    set(SpecUtils_ENABLE_D3_CHART ON)
    set(USE_SPECTRUM_CHART_D3 ON)
endif(BUILD_FOR_WEB_DEPLOYMENT)

if(BUILD_AS_LOCAL_SERVER)
    set(PERFORM_DEVELOPER_CHECKS ON)
endif(BUILD_AS_LOCAL_SERVER)

#check to make sure exactly one target is specified
set(NUM_TARGETS_SPECIFIED 0)

#if( BUILD_AS_OSX_APP )
#  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
#endif()

if(BUILD_AS_ELECTRON_APP)
    MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if(BUILD_FOR_WEB_DEPLOYMENT)
    MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if(BUILD_AS_LOCAL_SERVER)
    MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if(BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE)
    MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if(BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)
    MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if(ANDROID OR IOS)
    MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

if(NOT ${NUM_TARGETS_SPECIFIED} EQUAL "1")
    message(
        FATAL
        "You must specify exactly one target, you specified ${NUM_TARGETS_SPECIFIED}"
    )
endif()

#try to link to static libraries whenever possible
#  doing this to boost libraries only adds about 0.6 megabytes
if(TRY_TO_STATIC_LINK)
    set(Boost_USE_STATIC_LIBS ON)
    if(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif(WIN32)
else(TRY_TO_STATIC_LINK)
    set(Boost_USE_STATIC_LIBS OFF)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif(TRY_TO_STATIC_LINK)

INCLUDE(cmake/FindWt.cmake)
INCLUDE(cmake/WtFindZlib.txt) #ToDo: switch to using CMake privided FindZLIB.cmake, e.g. use find_package( ZLIB REQUIRED )

SET( Boost_ADDITIONAL_VERSIONS "1.55" "1.57" "1.65" "1.65_1" )
#find_package( Boost 1.41.0 COMPONENTS thread date_time system filesystem program_options regex random chrono atomic )

INCLUDE(cmake/WtFindBoost-cmake.txt)
if(NOT BOOST_WT_FOUND)
    SET(ERR
        "Could not find a boost installation in "   ${BOOST_PREFIX}   ".\n\n"
        "There are two methods in Wt to find boost:"
    )
    MESSAGE(FATAL_ERROR ${ERR})
endif(NOT BOOST_WT_FOUND)

set(sources
    src/InterSpecApp.cpp
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
    src/InterSpecUser.cpp
    src/CompactFileManager.cpp
    src/DetectorEdit.cpp
    src/OneOverR2Calc.cpp
    src/UnitsConverterTool.cpp
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
    src/FluxTool.cpp
    src/NativeFloatSpinBox.cpp
    src/EnergyCalTool.cpp
    src/EnergyCal.cpp
    src/EnergyCalGraphical.cpp
    src/EnergyCalAddActions.cpp
    src/EnergyCalPreserveWindow.cpp
    src/EnergyCalMultiFile.cpp
    src/SimpleDialog.cpp
    src/ShowRiidInstrumentsAna.cpp
    js/CanvasForDragging.js
    js/SpectrumChart.js
    js/InterSpec.js
)

set(headers
    InterSpec/InterSpec_config.h.in
    InterSpec/InterSpecApp.h
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
    InterSpec/InterSpecUser.h
    InterSpec/CompactFileManager.h
    InterSpec/DetectorEdit.h
    InterSpec/OneOverR2Calc.h
    InterSpec/UnitsConverterTool.h
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
    InterSpec/FluxTool.h
    InterSpec/NativeFloatSpinBox.h
    InterSpec/EnergyCalTool.h
    InterSpec/EnergyCal.h
    InterSpec/EnergyCalGraphical.h
    InterSpec/EnergyCalAddActions.h
    InterSpec/EnergyCalPreserveWindow.h
    InterSpec/EnergyCalMultiFile.h
    InterSpec/SimpleDialog.h
    InterSpec/ShowRiidInstrumentsAna.h
)

if(USE_DB_TO_STORE_SPECTRA)
    list(APPEND sources src/DbStateBrowser.cpp src/DbFileBrowser.cpp)
    list(APPEND headers InterSpec/DbStateBrowser.h InterSpec/DbFileBrowser.h)
endif(USE_DB_TO_STORE_SPECTRA)

if(ALLOW_URL_TO_FILESYSTEM_MAP)
    list(APPEND sources src/DbToFilesystemLink.cpp)
    list(APPEND headers InterSpec/DbToFilesystemLink.h)
endif(ALLOW_URL_TO_FILESYSTEM_MAP)

if(NOT ANDROID AND NOT IOS)
    list(APPEND sources src/FileDragUploadResource.cpp)
    list(APPEND headers InterSpec/FileDragUploadResource.h)
endif(NOT ANDROID AND NOT IOS)

if(INCLUDE_ANALYSIS_TEST_SUITE)
    list(APPEND sources src/SpectrumViewerTester.cpp)
    list(APPEND headers InterSpec/SpectrumViewerTester.h)
endif(INCLUDE_ANALYSIS_TEST_SUITE)

if(DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK)
    list(APPEND SRCS src/ChartToImageResource.cpp)
    list(APPEND HEADRS InterSpec/ChartToImageResource.h)
endif(DECAY_CHART_ADD_IMAGE_DOWNLOAD_LINK)

if(BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)
    list(APPEND sources testing/developcode.cpp src/DbUserState.cpp)
    list(APPEND headers testing/developcode.h InterSpec/DbUserState.h)
endif(BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)

if(USE_OSX_NATIVE_MENU)
    list(APPEND sources target/osx/NativeMenu.mm)
    list(APPEND headers target/osx/NativeMenu.h)
endif(USE_OSX_NATIVE_MENU)

if(IOS)
    list(APPEND sources target/ios/InterSpec/FileHandling.mm)
    list(APPEND headers target/ios/InterSpec/FileHandling.h)
endif(IOS)

if(USE_GOOGLE_MAP)
    list(APPEND sources src/GoogleMap.cpp)
    list(APPEND headers InterSpec/GoogleMap.h)
endif(USE_GOOGLE_MAP)

if(USE_SEARCH_MODE_3D_CHART)
    list(APPEND sources src/SearchMode3DChart.cpp src/SearchMode3DDataModel.cpp)
    list(
        APPEND
        headers
        InterSpec/SearchMode3DChart.h
        InterSpec/SearchMode3DDataModel.h
    )
endif(USE_SEARCH_MODE_3D_CHART)

if(USE_FEATURE_MARKER_WIDGET)
    list(APPEND sources src/FeatureMarkerWidget.cpp)
    list(APPEND headers InterSpec/FeatureMarkerWidget.h)
endif(USE_FEATURE_MARKER_WIDGET)

if(USE_TERMINAL_WIDGET)
    list(
        APPEND
        sources
        src/TerminalModel.cpp
        src/TerminalWidget.cpp
        js/TerminalWidget.js
    )
    list(APPEND headers InterSpec/TerminalModel.h InterSpec/TerminalWidget.h)
endif(USE_TERMINAL_WIDGET)

if(SUPPORT_ZIPPED_SPECTRUM_FILES)
    list(APPEND sources src/ZipArchive.cpp)
    list(APPEND headers InterSpec/ZipArchive.h)
endif(SUPPORT_ZIPPED_SPECTRUM_FILES)

if(USE_SIMPLE_NUCLIDE_ASSIST)
    list(APPEND sources src/SimpleNuclideAssist.cpp)
    list(APPEND headers InterSpec/SimpleNuclideAssist.h)
endif(USE_SIMPLE_NUCLIDE_ASSIST)

if(USE_SPECRUM_FILE_QUERY_WIDGET)
    if(
        BUILD_AS_OSX_APP
        OR BUILD_AS_ELECTRON_APP
        OR BUILD_AS_LOCAL_SERVER
        OR BUILD_AS_UNIT_TEST_SUITE
    )
        list(
            APPEND
            sources
            src/SpecFileQuery.cpp
            src/SpecFileQueryWidget.cpp
            src/SpecFileQueryDbCache.cpp
        )
        list(
            APPEND
            headers
            InterSpec/SpecFileQuery.h
            InterSpec/SpecFileQueryWidget.h
            InterSpec/SpecFileQueryDbCache.h
        )
        list(APPEND sources js/SpecFileQueryWidget.js)
    else()
        message(
            FATAL_ERROR
            "The spectrum file query widget can not be enabled for web deployements."
        )
    endif()
endif(USE_SPECRUM_FILE_QUERY_WIDGET)

#Copy D3 resources into InterSpec_resources directory at compile time
include(cmake/DeployJsAndCss.cmake)

list(APPEND sources external_libs/SpecUtils/d3_resources/d3.v3.min.js)
deploy_js_resource("${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/d3_resources/d3.v3.min.js"
                    "${CMAKE_CURRENT_BINARY_DIR}/InterSpec_resources/d3.v3.min.js"
)

if(USE_SPECTRUM_CHART_D3)
    if(NOT RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE)
        MESSAGE(FATAL_ERROR     "You must enable RENDER_REFERENCE_PHOTOPEAKS_SERVERSIDE when enabling USE_SPECTRUM_CHART_D3")
    endif()

    list(
        APPEND
        sources
        src/D3SpectrumDisplayDiv.cpp
        external_libs/SpecUtils/d3_resources/SpectrumChartD3.js
        external_libs/SpecUtils/d3_resources/SpectrumChartD3.css
        external_libs/SpecUtils/d3_resources/SpectrumChartD3StandAlone.css
        src/D3TimeChart.cpp
    )
    list( APPEND headers InterSpec/D3SpectrumDisplayDiv.h InterSpec/D3TimeChart.h )

    list( APPEND OTHER_SUPPORT_FILES ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/D3TimeChart.js )

    deploy_js_resource("${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/d3_resources/SpectrumChartD3.js"
         "${CMAKE_CURRENT_BINARY_DIR}/InterSpec_resources/SpectrumChartD3.js"
    )

    set(CSS_RESOURCES_TO_COPY
        "SpectrumChartD3.css"
        "SpectrumChartD3StandAlone.css"
    )
    foreach(_file ${CSS_RESOURCES_TO_COPY})
        deploy_css_resource("${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/d3_resources/${_file}"
          "${CMAKE_CURRENT_BINARY_DIR}/InterSpec_resources/${_file}"
        )
    endforeach()
endif(USE_SPECTRUM_CHART_D3)

find_path(
    WT_RESOURCES_DIRECTORY
    form.css
    ${Wt_INCLUDE_DIR}/../share/Wt/resources
)

#Set SpecUtils compile options
set(SpecUtils_REBIN_FILES_TO_SINGLE_BINNING OFF
    CACHE BOOL "Set Rebin spectrum to single energy cal"
)
set(SpecUtils_NO_BOOST_LIB OFF
    CACHE BOOL "Use boost for threading primitives in SpecUtils"
)
set(SpecUtils_USE_WT_THREADPOOL ON CACHE BOOL "")
set(SpecUtils_EXTERNALLY_DEFINED_LOG_MESSAGE ON CACHE BOOL "")
set(SpecUtils_BUILD_UNIT_TESTS OFF CACHE BOOL "")
set(SpecUtils_BUILD_REGRESSION_TEST OFF CACHE BOOL "")
set(SpecUtils_PYTHON_BINDINGS OFF CACHE BOOL "")
set(SpecUtils_JAVA_SWIG OFF CACHE BOOL "")
set(SpecUtils_D3_SUPPORT_FILE_STATIC OFF CACHE BOOL "")
set(SpecUtils_D3_SCRIPTS_RUNTIME_DIR "InterSpec_resources" CACHE STRING "")

#add in all the subdirectories the project will always use
add_subdirectory(external_libs/Minuit2)
add_subdirectory(external_libs/Cuba-3.0)
add_subdirectory(external_libs/SpecUtils)
add_subdirectory(external_libs/SandiaDecay)
add_subdirectory(external_libs/muparserx-4.0.7)

#We should have figured out all the header and source files we need by here, so
#  lets add the executable to the build
if(IOS)
    list(APPEND sources src/InterSpecServer.cpp)
    list(APPEND headers InterSpec/InterSpecServer.h)
    add_library(InterSpecLib STATIC ${sources} ${headers})
elseif(ANDROID)
    list(APPEND headers target/android/AndroidUtils.hpp)
    add_library(
        InterSpecLib
        SHARED
        main.cpp
        ${sources}
        ${headers}
    )
elseif(BUILD_AS_OSX_APP)
    list(APPEND sources src/InterSpecServer.cpp target/osx/macOsUtils.mm)
    list(APPEND headers InterSpec/InterSpecServer.h target/osx/macOsUtils.h)

    SET(MACOSX_BUNDLE_INFO_STRING   "InterSpec - Sandia National Lab, Will Johnson")
    SET(MACOSX_BUNDLE_ICON_FILE InterSpec )
    SET(MACOSX_BUNDLE_LONG_VERSION_STRING "" )
    SET(MACOSX_BUNDLE_BUNDLE_NAME "InterSpec" )
    SET(MACOSX_BUNDLE_SHORT_VERSION_STRING "InterSpec" )
    SET(MACOSX_BUNDLE_BUNDLE_VERSION "0" )
    SET(MACOSX_BUNDLE_COPYRIGHT "Sandia National Labs, Will Johnson" )
    SET(MACOSX_BUNDLE_GUI_IDENTIFIER "gov.sandia.macOS.InterSpec" )

    SET( ${PRODUCT_NAME} "InterSpec" )
    SET( MACOSX_BUNDLE_BUNDLE_NAME "InterSpec" )

    FIND_LIBRARY(WEBKIT_LIBRARY WebKit)
    FIND_LIBRARY(COCOA_LIBRARY Cocoa)
    FIND_LIBRARY(COREDATA_LIBRARY CoreData)
    FIND_LIBRARY(APPKIT_LIBRARY AppKit)
    FIND_LIBRARY(FOUNDATION_LIBRARY Foundation)

    MARK_AS_ADVANCED(WEBKIT_LIBRARY COCOA_LIBRARY COREDATA_LIBRARY APPKIT_LIBRARY FOUNDATION_LIBRARY)
    add_library(InterSpecLib STATIC ${sources} ${headers})

    target_link_libraries(
        InterSpecLib
        PUBLIC
            ${WEBKIT_LIBRARY}
            ${COCOA_LIBRARY}
            ${APPKIT_LIBRARY}
            ${COREDATA_LIBRARY}
            ${FOUNDATION_LIBRARY}
    )

    set(app_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/target/osx/InterSpec.icns
        "target/osx/en.lproj/MainMenu.xib"
        "target/osx/AppDelegate.mm"
        "target/osx/main.m"
    )

    set(app_headers ${headers} "target/osx/AppDelegate.h")

    # Xcode doesnt like having two bundles for some reason, so we will add the
    # QuickLook project as an external project and add it to the final application
    # at InterSpec.app/Contents/Library/QuickLook
    include(ExternalProject)
    ExternalProject_Add(SpecFilePreview
          BINARY_DIR   "${CMAKE_BINARY_DIR}/SpecFilePreview"
          URL ${CMAKE_CURRENT_SOURCE_DIR}/target/macOsQuickLook/SpecFilePreview/SpecFilePreview/
          CMAKE_GENERATOR   "Unix Makefiles"
          CMAKE_CACHE_ARGS -DWt_INCLUDE_DIR:PATH=${Wt_INCLUDE_DIR}
                           -DBOOST_ROOT:PATH=${Boost_INCLUDE_DIR}/..
                           -DBOOST_INCLUDEDIR:PATH=${Boost_INCLUDE_DIR}
                           -DCMAKE_BUILD_TYPE:STRING=RELEASE
                           -DSpecUtils_ENABLE_D3_CHART:BOOL=OFF
                           INSTALL_COMMAND   ""
    )

    # Add dummy command to let cmake know SpecFilePreview.qlgenerator isnt available
    # at configure time so it wont complain
    add_custom_command(
        OUTPUT
            "${CMAKE_CURRENT_BINARY_DIR}/SpecFilePreview/SpecFilePreview.qlgenerator"
        COMMAND ""
    )

    #In order to copy all the needed resources to the "Archive" builds of Xcode,
    #  we need to specify each individual file, and later set their MACOSX_PACKAGE_LOCATION properties...
    file(GLOB_RECURSE DATA_RES_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/data/*")
    file(GLOB_RECURSE WT_RES_SOURCES "${WT_RESOURCES_DIRECTORY}/*")
    file(
        GLOB_RECURSE INTERSPEC_RES_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/*"
    )
    file(
        GLOB_RECURSE EXAMPLE_RES_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/example_spectra/*"
    )
    SET(EN_LOCALIZATION_FILES   "${CMAKE_CURRENT_SOURCE_DIR}/target/osx/en.lproj/InfoPlist.strings"   "${CMAKE_CURRENT_SOURCE_DIR}/target/osx/en.lproj/Credits.rtf")

    find_program(IBTOOL ibtool HINTS "/usr/bin" "${OSX_DEVELOPER_ROOT}/usr/bin")
    if(${IBTOOL} STREQUAL "IBTOOL-NOTFOUND")
        message(SEND_ERROR "could fine ibtool to compile the .xib files")
    endif()

    add_custom_command(
        OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/MainMenu.nib"
        COMMAND
            ${IBTOOL} --errors --warnings --notices --output-format
            human-readable-text --compile
            "${CMAKE_CURRENT_BINARY_DIR}/MainMenu.nib"
            "${CMAKE_CURRENT_SOURCE_DIR}/target/osx/en.lproj/MainMenu.xib"
        COMMENT "Compiling MainMenu.xib"
    )

    ADD_EXECUTABLE(InterSpecExe MACOSX_BUNDLE ${app_sources} ${app_headers}
                      ${DATA_RES_SOURCES} ${WT_RES_SOURCES} ${INTERSPEC_RES_SOURCES} ${EXAMPLE_RES_SOURCES} ${EN_LOCALIZATION_FILES}
                      "${CMAKE_CURRENT_BINARY_DIR}/MainMenu.nib"
                      "${CMAKE_CURRENT_BINARY_DIR}/SpecFilePreview/SpecFilePreview.qlgenerator"
    )
    SET_TARGET_PROPERTIES( InterSpecExe PROPERTIES MACOSX_BUNDLE TRUE )
    #set_target_properties( InterSpecExe PROPERTIES SUFFIX "")
    set_target_properties(InterSpecExe PROPERTIES OUTPUT_NAME "InterSpec")

    target_link_libraries(InterSpecExe PUBLIC InterSpecLib)

    SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_SOURCE_DIR}/target/osx/InterSpec.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)

    add_dependencies(InterSpecExe SpecFilePreview)

    #Set resources MACOSX_PACKAGE_LOCATION property so it is where we want in the app bundle
    #ToDo: combine 5 loops below into a single loop...
    # ToDo: Try to use: install(DIRECTORY InterSpec_resources DESTINATION $<TARGET_FILE_DIR:InterSpecExe>/../Resources/ )
    # ToDo: Try using (something like):
    #                 add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory
    #                                    "${CMAKE_SOURCE_DIR}/example_spectra" "$<TARGET_FILE_DIR:InterSpecExe>/Contents/Resources/")

    # set( RES_SOURCES_DIRS DATA_RES_SOURCES WT_RES_SOURCES INTERSPEC_RES_SOURCES EXAMPLE_RES_SOURCES EN_LOCALIZATION_FILES )
    # set( RES_DEST_DIRS data resources InterSpec_resources example_spectra en.lproj )

    set_property(
        SOURCE
        ${CMAKE_CURRENT_BINARY_DIR}/SpecFilePreview/SpecFilePreview.qlgenerator
        PROPERTY
        MACOSX_PACKAGE_LOCATION
        "Library/QuickLook/"
    )

    foreach(RES_FILE ${DATA_RES_SOURCES})
        #Get the relative path from the data-folder to the particular file
        file(
            RELATIVE_PATH
            RES_PATH
            "${CMAKE_CURRENT_SOURCE_DIR}/data"
            ${RES_FILE}
        )
        get_filename_component(RES_PATH ${RES_PATH} DIRECTORY)
        #Set it's location inside the app package (under Resources)
        set_property(
            SOURCE
            ${RES_FILE}
            PROPERTY
            MACOSX_PACKAGE_LOCATION
            "Resources/data/${RES_PATH}"
        )
    endforeach(RES_FILE)

    foreach(RES_FILE ${WT_RES_SOURCES})
        file(RELATIVE_PATH RES_PATH "${WT_RESOURCES_DIRECTORY}" ${RES_FILE})
        get_filename_component(RES_PATH ${RES_PATH} DIRECTORY)
        set_property(
            SOURCE
            ${RES_FILE}
            PROPERTY
            MACOSX_PACKAGE_LOCATION
            "Resources/resources/${RES_PATH}"
        )
    endforeach(RES_FILE)

    foreach(RES_FILE ${INTERSPEC_RES_SOURCES})
        file(
            RELATIVE_PATH
            RES_PATH
            "${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources"
            ${RES_FILE}
        )
        get_filename_component(RES_PATH ${RES_PATH} DIRECTORY)
        set_property(
            SOURCE
            ${RES_FILE}
            PROPERTY
            MACOSX_PACKAGE_LOCATION
            "Resources/InterSpec_resources/${RES_PATH}"
        )
    endforeach(RES_FILE)

    foreach(RES_FILE ${EXAMPLE_RES_SOURCES})
        file(
            RELATIVE_PATH
            RES_PATH
            "${CMAKE_CURRENT_SOURCE_DIR}/example_spectra"
            ${RES_FILE}
        )
        get_filename_component(RES_PATH ${RES_PATH} DIRECTORY)
        set_property(
            SOURCE
            ${RES_FILE}
            PROPERTY
            MACOSX_PACKAGE_LOCATION
            "Resources/example_spectra/${RES_PATH}"
        )
    endforeach(RES_FILE)

    foreach(RES_FILE ${EN_LOCALIZATION_FILES})
        set_property(
            SOURCE
            ${RES_FILE}
            PROPERTY
            MACOSX_PACKAGE_LOCATION
            "Resources/en.lproj/"
        )
    endforeach(RES_FILE)

    set_property(
        SOURCE
        "${CMAKE_CURRENT_BINARY_DIR}/MainMenu.nib"
        PROPERTY
        MACOSX_PACKAGE_LOCATION
        "Resources/en.lproj/"
    )

    #set MACOS_BUNDLE_VERSION_NUMBER for Info.plist.template to populate cooresponding value in Xcode
    #math(EXPR MACOS_BUNDLE_VERSION_NUMBER "100*(${PROJECT_VERSION_MAJOR}-1) + 10*${PROJECT_VERSION_MINOR} + ${PROJECT_VERSION_PATCH}" )
    set(MACOS_BUNDLE_VERSION_NUMBER 9)

    set_target_properties(
        InterSpecExe
        PROPERTIES
        MACOSX_BUNDLE_INFO_PLIST
        "${CMAKE_CURRENT_SOURCE_DIR}/target/osx/Info.plist.template"
    )

    set_target_properties(
        InterSpecExe
        PROPERTIES
        XCODE_ATTRIBUTE_CODE_SIGN_ENTITLEMENTS
        "${CMAKE_CURRENT_SOURCE_DIR}/target/osx/InterSpec.entitlements"
    )

    # With out this next line, Xcode will (silently) fail to build Archive builds
    set_target_properties(
        InterSpecExe
        PROPERTIES
        XCODE_ATTRIBUTE_INSTALL_PATH
        "/Applications"
    )

    #Figure out the development team ID, and set this info for Xcode so it will default to signing the app
    #  This is probably a little brittle on other peoples machined
    execute_process(
        COMMAND
            /usr/bin/security find-identity -v -p codesigning OUTPUT_VARIABLE
            dev_team_ids
    )
    STRING(REGEX REPLACE ";" "\\\\;" dev_team_ids "${dev_team_ids}")
    STRING(REGEX REPLACE "\n" ";" dev_team_ids "${dev_team_ids}")
    LIST(FILTER dev_team_ids INCLUDE REGEX "Developer ID Application:.+" )
    STRING(REGEX MATCH "\\(.+\\)" mac_team_id ${dev_team_ids} )
    STRING(REGEX REPLACE "\\(|\\)" "" mac_team_id "${mac_team_id}")
    message("Using Mac Developer Team ID='${mac_team_id}'")

    set_target_properties(
        InterSpecExe
        PROPERTIES
        XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY
        "Mac Developer"
    )
    set_target_properties(
        InterSpecExe
        PROPERTIES
        XCODE_ATTRIBUTE_DEVELOPMENT_TEAM
        "${mac_team_id}"
    )
    set_target_properties(
        InterSpecExe
        PROPERTIES
        XCODE_ATTRIBUTE_ENABLE_HARDENED_RUNTIME
        "YES"
    )
elseif(BUILD_AS_UNIT_TEST_SUITE)
    add_library(InterSpecLib STATIC ${sources} ${headers})
    add_executable(InterSpecExe ${GUI_TYPE} main.cpp)
    target_link_libraries(InterSpecExe PUBLIC InterSpecLib)
    set_target_properties(InterSpecExe PROPERTIES OUTPUT_NAME "InterSpec")
elseif(BUILD_AS_ELECTRON_APP)
    list(
        APPEND
        sources
        src/InterSpecServer.cpp
        target/electron/ElectronUtils.cpp
    )
    
    list(
        APPEND
        headers
        InterSpec/InterSpecServer.h
        target/electron/ElectronUtils.h
    )
    
    if( NOT USE_ELECTRON_NATIVE_MENU )
      list( APPEND sources js/AppHtmlMenu.js )
    endif( NOT USE_ELECTRON_NATIVE_MENU )
    
    if( INTERSPEC_LIBRARY_STATIC )
        add_library(InterSpecLib STATIC ${sources} ${headers})
    else( INTERSPEC_LIBRARY_STATIC )
        add_library(InterSpecLib SHARED ${sources} ${headers})
    endif( INTERSPEC_LIBRARY_STATIC )

    #set_property( TARGET InterSpecLib PROPERTY POSITION_INDEPENDENT_CODE ON )
else()
    add_library(InterSpecLib STATIC ${sources} ${headers})
    add_executable(InterSpecExe ${GUI_TYPE} main.cpp)
    target_link_libraries(InterSpecExe PUBLIC InterSpecLib)
    #set_target_properties( InterSpecExe PROPERTIES SUFFIX "")
    set_target_properties(InterSpecExe PROPERTIES OUTPUT_NAME "InterSpec")

#    add_executable(peakFitTimingDev testing/peakFitTimingDev.cpp)
#    target_link_libraries(peakFitTimingDev PUBLIC InterSpecLib)
endif(IOS)

set_target_properties(InterSpecLib PROPERTIES PREFIX "")
set_target_properties(InterSpecLib PROPERTIES OUTPUT_NAME "InterSpec")

target_link_libraries(
    InterSpecLib
    PUBLIC
    SpecUtils
    SandiaDecay
    minuit2
    Cuba-3.0
    muparserx
)

if(USE_SPECRUM_FILE_QUERY_WIDGET OR BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT)
    add_subdirectory(external_libs/pugixml-1.9)
    target_link_libraries(InterSpecLib PUBLIC pugixml)
endif()

if(
    BUILD_AS_OFFLINE_ANALYSIS_TEST_SUITE
    OR
    BUILD_AS_UNIT_TEST_SUITE
    OR
    BUILD_AS_COMMAND_LINE_CODE_DEVELOPMENT
)
    target_link_libraries(InterSpecLib PUBLIC ${Wt_TEST_LIBRARY})
endif()

if(USE_MYSQL_DB EQUAL USE_SQLITE3_DB)
    MESSAGE(FATAL_ERROR   "Exactly one type of database (MySQL or SQLITE3) must be chosen")
endif()

set(RAPID_XML_INC_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/3rdparty
)

target_include_directories(
    InterSpecLib
    PUBLIC
    ${Wt_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
    ${RAPID_XML_INC_DIR}
    external_libs
)

if(IOS)
elseif(WIN32)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_link_libraries(InterSpecLib PUBLIC ${Wt_DEBUG_LIBRARIES})
    else(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_link_libraries(InterSpecLib PUBLIC ${Wt_LIBRARIES})
    endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

    target_link_directories(InterSpecLib PUBLIC ${Boost_LIBRARY_DIR})

    if(Boost_CHRONO_FOUND)
        target_link_libraries(InterSpecLib PUBLIC ${Boost_CHRONO_LIBRARY})
    endif(Boost_CHRONO_FOUND)

    if(Boost_TIMER_FOUND)
        target_link_libraries(InterSpecLib PUBLIC ${Boost_TIMER_LIBRARY})
    endif(Boost_TIMER_FOUND)
elseif(BUILD_AS_ELECTRON_APP)
    #Not sure why Wt_LIBRARIES isnt filled out here, but its not, so will manually specify Wt libraries for the moment
    target_link_libraries(
        InterSpecLib
        PUBLIC
        ${Wt_LIBRARY}
        ${Wt_HTTP_LIBRARY}
        ${Wt_DBO_LIBRARY}
        ${Wt_DBOSQLITE3_LIBRARY}
    )

    if(INCLUDE_ANALYSIS_TEST_SUITE)
        target_link_libraries(InterSpecLib PUBLIC ${Wt_TEST_LIBRARY})
    endif(INCLUDE_ANALYSIS_TEST_SUITE)
elseif(YES)
    target_link_libraries(InterSpecLib PUBLIC ${Wt_LIBRARIES})
endif(IOS)

if(NOT IOS)
    target_link_libraries(
        InterSpecLib
        PUBLIC
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_REGEX_LIBRARY}
        ${Boost_THREAD_LIBRARY}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_RANDOM_LIBRARY}
        ${Boost_PROGRAM_OPTIONS_LIBRARY}
        ${Boost_DATE_TIME_LIBRARY}
    )

    find_package(Threads REQUIRED)
    target_link_libraries(InterSpecLib PUBLIC Threads::Threads)
endif(NOT IOS)

if(APPLE AND NOT IOS)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.12")
endif()

if( BUILD_AS_LOCAL_SERVER AND (CMAKE_GENERATOR STREQUAL "Xcode") )
    
    # Setting Xcode working directory requires cmake 3.17
    if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17.0" )
        set_target_properties(InterSpecExe PROPERTIES
            XCODE_GENERATE_SCHEME TRUE
            XCODE_SCHEME_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        )
    endif( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17.0" )
    
    # Xcode stopped generating debug sybols by default recently, so force enable,
    if( CMAKE_BUILD_TYPE MATCHES Debug )
      target_compile_options(SpecUtils PUBLIC $<$<CXX_COMPILER_ID:Clang>:-g> )
      target_compile_options(SandiaDecay PUBLIC $<$<CXX_COMPILER_ID:Clang>:-g> )
      target_compile_options(minuit2 PUBLIC $<$<CXX_COMPILER_ID:Clang>:-g> )
      target_compile_options(InterSpecLib PUBLIC $<$<CXX_COMPILER_ID:Clang>:-g> )
      #Cuba-3.0, muparserx
    endif( CMAKE_BUILD_TYPE MATCHES Debug )
endif( BUILD_AS_LOCAL_SERVER AND (CMAKE_GENERATOR STREQUAL "Xcode") )

if( BUILD_AS_LOCAL_SERVER AND (CMAKE_GENERATOR MATCHES "Visual Studio") )
    #set_target_properties(InterSpecExe PROPERTIES
    #    VS_DEBUGGER_COMMAND "debug_command"
    #    VS_DEBUGGER_COMMAND_ARGUMENTS "debug_arguments")
endif( BUILD_AS_LOCAL_SERVER AND (CMAKE_GENERATOR MATCHES "Visual Studio") )


#20180104: on my mac I tried compiling Wt with support for PDF (libharu), but only
#  had the static libhpdfs, so I had to add bellow to link to this
FIND_LIBRARY( HPDF_LIBRARY NAME hpdfs HINTS ${Wt_INCLUDE_DIR}/../lib )
if(NOT HPDF_LIBRARY)
    message(
        "Could not find libharu - if Wt was built with PDF support, you will get link errors"
    )
else()
    target_link_libraries(InterSpecLib PUBLIC ${HPDF_LIBRARY})
    FIND_LIBRARY( PNG_LIBRARY NAME png HINTS ${Wt_INCLUDE_DIR}/../lib )
    if(NOT PNG_LIBRARY)
        message(FATAL "Could not find libpng - requred for libharu")
    else()
        target_link_libraries(InterSpecLib PUBLIC ${PNG_LIBRARY})
    endif()
endif()

#bellow is what I need to get things working 20131120 for statically linked
if(TRY_TO_STATIC_LINK)
    if(APPLE AND NOT IOS)
        FIND_LIBRARY(LIBZ_LIBRARY z)

        target_link_libraries(
            InterSpecLib
            PUBLIC
            ${Boost_RANDOM_LIBRARY}
            ${Boost_PROGRAM_OPTIONS_LIBRARY}
            ${LIBZ_LIBRARY}
            libSystem.B.dylib
            /System/Library/Frameworks/Carbon.framework/Versions/A/Carbon
            /System/Library/Frameworks/Foundation.framework/Versions/C/Foundation
        )

        if(USE_OSX_NATIVE_MENU)
            target_link_libraries(
                InterSpecLib
                PUBLIC
                /System/Library/Frameworks/AppKit.framework/Versions/C/AppKit
            )
        endif(USE_OSX_NATIVE_MENU)
    endif(APPLE AND NOT IOS)

    target_link_libraries(InterSpecLib PUBLIC ${CMAKE_DL_LIBS})

    if(USE_MYSQL_DB)
        include(cmake/WtFindMysql.txt)
        target_link_libraries(InterSpecLib PUBLIC ${MYSQL_LIBRARIES})
    endif(USE_MYSQL_DB)

    if(Wt_FCGI_LIBRARY)
        include(cmake/WtFindFcgi.txt)
        target_link_libraries(InterSpecLib PUBLIC ${FCGI_LIBRARIES})
    endif(Wt_FCGI_LIBRARY)

    include(cmake/WtFindSsl.txt)
    if(SSL_FOUND)
        target_link_libraries(InterSpecLib PUBLIC ${SSL_LIBRARIES})
    endif(SSL_FOUND)

    #link to zlib since Wt might
    if("${CMAKE_SYSTEM}" MATCHES "Linux" AND NOT SUPPORT_ZIPPED_SPECTRUM_FILES)
        target_link_libraries(InterSpecLib PUBLIC ${ZLIB_LIBRARIES})
    endif(
        "${CMAKE_SYSTEM}" MATCHES "Linux"
        AND
        NOT SUPPORT_ZIPPED_SPECTRUM_FILES
    )

    if("${CMAKE_SYSTEM}" MATCHES "Linux" AND NOT ANDROID)
        find_library(RT_LIB rt)
        if(RT_LIB EQUAL RT_LIB-NOTFOUND)
            message(FATAL "could not find rt lib")
        endif()

        target_link_libraries(InterSpecLib PUBLIC ${RT_LIB})

        #target_link_libraries( InterSpecExe -static-libgcc -static-libstdc++)
        #Next line statically links to pthread, copying whole archive into exe
        #set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-static -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    endif("${CMAKE_SYSTEM}" MATCHES "Linux" AND NOT ANDROID)
endif(TRY_TO_STATIC_LINK)

if(BUILD_AS_UNIT_TEST_SUITE)
    target_link_libraries(
        InterSpecLib
        PUBLIC
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    )
endif(BUILD_AS_UNIT_TEST_SUITE)

if(PERFORM_DEVELOPER_CHECKS)
    target_link_libraries(InterSpecLib PUBLIC ${Boost_DATE_TIME_LIBRARY})
endif(PERFORM_DEVELOPER_CHECKS)

include(CheckIncludeFileCXX)
set(CMAKE_REQUIRED_INCLUDES ${Wt_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )

if(Wt_DBOSQLITE3_LIBRARY)
    # check_include_file_cxx(...) doesnt look to work w/ apple clang 12...
    #check_include_file_cxx("Wt/Dbo/backend/Sqlite3" HAS_WT_DBO_SQLITE3)
    find_file( HAS_WT_DBO_SQLITE3 "Sqlite3" PATHS "${Wt_INCLUDE_DIR}/Wt/Dbo/backend/" NO_DEFAULT_PATH )
endif(Wt_DBOSQLITE3_LIBRARY)

if(Wt_MYSQL_LIBRARY)
    #check_include_file_cxx("${Wt_INCLUDE_DIR}/Wt/Dbo/backend/MySQL" HAS_WT_DBO_MYSQL)
    find_file( HAS_WT_DBO_MYSQL "MySQL" PATHS "${Wt_INCLUDE_DIR}/Wt/Dbo/backend/" NO_DEFAULT_PATH )
endif(Wt_MYSQL_LIBRARY)

INCLUDE( CheckIncludeFiles )

if(ZLIB_FOUND AND Boost_IOSTREAMS_FOUND)
    #CHECK_INCLUDE_FILES( "zlib.h" HAS_ZLIB_SUPPORT )
    target_link_libraries(
        InterSpecLib
        PUBLIC
        ${ZLIB_LIBRARIES}
        ${Boost_IOSTREAMS_LIBRARY}
    )
endif(ZLIB_FOUND AND Boost_IOSTREAMS_FOUND)

if(SUPPORT_ZIPPED_SPECTRUM_FILES)
    #CMake seems to fail to find zlib on ios/android (even using find_package(ZLIB REQUIRED)), but we know it will be in
    #  the systems path anyway, so just blindly link to zlib
    if(NOT (IOS OR ANDROID))
        if(NOT ZLIB_FOUND)
            MESSAGE(FATAL_ERROR       "ZLIB required to enable SUPPORT_ZIPPED_SPECTRUM_FILES, but it wasnt found")
        else(NOT ZLIB_FOUND)
            target_include_directories(InterSpecLib PUBLIC ${ZLIB_INCLUDE_DIRS})
            target_link_libraries(InterSpecLib PUBLIC ${ZLIB_LIBRARIES})
        endif(NOT ZLIB_FOUND)
    else(NOT (IOS OR ANDROID))
        target_link_libraries(InterSpecLib PUBLIC z)
    endif(NOT (IOS OR ANDROID))
endif(SUPPORT_ZIPPED_SPECTRUM_FILES)

if(Wt_MYSQL_LIBRARY AND BUILD_FOR_WEB_DEPLOYMENT AND USE_MYSQL_DB)
    set(DATABASE_PASSWORD_FILE "/path/to/passwords/databases.xml"
        CACHE STRING
        "Location of the XML passwords file for connecting to MySQL database"
    )
endif(Wt_MYSQL_LIBRARY AND BUILD_FOR_WEB_DEPLOYMENT AND USE_MYSQL_DB)

if(USE_MYSQL_DB AND NOT HAS_WT_DBO_MYSQL)
    MESSAGE(FATAL_ERROR   "You do not have the MySQL Wt Dbo backend installed, but you have asked to use it")
endif(USE_MYSQL_DB AND NOT HAS_WT_DBO_MYSQL)

if(USE_SQLITE3_DB AND NOT HAS_WT_DBO_SQLITE3 AND NOT WIN32)
    MESSAGE(FATAL_ERROR   "You do not have the SQLITE3 Wt Dbo backend installed, but you have asked to use it")
endif(USE_SQLITE3_DB AND NOT HAS_WT_DBO_SQLITE3 AND NOT WIN32)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec/InterSpec_config.h.in
    ${CMAKE_BINARY_DIR}/InterSpec_config.h
)

if(ANDROID)
    target_link_libraries(InterSpecLib PUBLIC log)
    target_include_directories(InterSpecLib PUBLIC target)
endif(ANDROID)

if(BUILD_AS_UNIT_TEST_SUITE)
    add_subdirectory(testing)
endif(BUILD_AS_UNIT_TEST_SUITE)

list( APPEND OTHER_SUPPORT_FILES
     ${CMAKE_CURRENT_SOURCE_DIR}/README.md
     ${CMAKE_CURRENT_SOURCE_DIR}/data/default_preferences.xml
     ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/themes/dark/dark.css
     ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/copyright_and_about.xml
     ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/use_instructions.xml
  )

FILE(GLOB CSS_FILES ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/*.css )
list(FILTER CSS_FILES EXCLUDE REGEX ".*SpectrumChartD3.*css")
list(APPEND OTHER_SUPPORT_FILES ${CSS_FILES})

list(
    APPEND
    OTHER_SUPPORT_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/DecayChainChart.js
)

FILE(GLOB HELP_FILES ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/*help.xml)
list(
    APPEND
    OTHER_SUPPORT_FILES
    ${HELP_FILES}
    ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/help.json
)

add_custom_target(interspec_resources SOURCES ${OTHER_SUPPORT_FILES})

if(BUILD_AS_ELECTRON_APP)
    set(ELECTRON_SUPPORT_FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/main.js
        ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/loading.html
        ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/README.md
        ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/package.json
        ${CMAKE_CURRENT_SOURCE_DIR}/target/electron/app/launch_options.json
    )

    add_custom_target(electron_resources SOURCES ${ELECTRON_SUPPORT_FILES})
endif(BUILD_AS_ELECTRON_APP)

#Copy files to the install directory
set(OTHER_FILES runLocalOn8080.sh)

foreach(_file ${OTHER_FILES})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_file}
        ${PROJECT_BINARY_DIR}/${_file}
        COPYONLY
    )
endforeach()

if(IOS)
    MESSAGE("CMAKE_CURRENT_SOURCE_DIR is ${CMAKE_CURRENT_SOURCE_DIR}")

    MESSAGE("Looking for form.css in ${CMAKE_CURRENT_SOURCE_DIR}/wt-3.3.1/resources")
    if(NOT WT_RESOURCES_DIRECTORY)
        FIND_PATH(WT_RESOURCES_DIRECTORY form.css HINTS ${CMAKE_CURRENT_SOURCE_DIR}/wt-3.3.1/resources NO_DEFAULT_PATHS NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_FIND_ROOT_PATH)
    endif(NOT WT_RESOURCES_DIRECTORY)
    MESSAGE("***** WT_RESOURCES_DIRECTORY is ${WT_RESOURCES_DIRECTORY}")

    MESSAGE("WT_RESOURCES_DIRECTORY IS ${WT_RESOURCES_DIRECTORY}")
endif(IOS)

if(NOT WT_RESOURCES_DIRECTORY)
    MESSAGE("Could not find the Wt resources path to copy into current directory")
else(NOT WT_RESOURCES_DIRECTORY)
    MESSAGE("Using Wt resources from ${WT_RESOURCES_DIRECTORY}")
    file(COPY ${WT_RESOURCES_DIRECTORY} DESTINATION ${PROJECT_BINARY_DIR}/)
endif(NOT WT_RESOURCES_DIRECTORY)

if(NOT ANDROID)
    SET( SUPPORT_DIRECTORIES InterSpec_resources data example_spectra )
endif(NOT ANDROID)

foreach(_dir ${SUPPORT_DIRECTORIES})
    if(WIN32)
        #        file(COPY ${PROJECT_SOURCE_DIR}/${_dir} DESTINATION ${PROJECT_BINARY_DIR}/ )
    else(WIN32)
        execute_process(
            COMMAND
                ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/${_dir}
                ${PROJECT_BINARY_DIR}/${_dir}
        )
    endif(WIN32)
endforeach()

if(INCLUDE_ANALYSIS_TEST_SUITE)
    if(WIN32)
        #        file(COPY ${PROJECT_SOURCE_DIR}/testing/analysis_tests DESTINATION ${PROJECT_BINARY_DIR}/analysis_tests )
    else(WIN32)
        execute_process(
            COMMAND
                ${CMAKE_COMMAND} -E create_symlink
                ${PROJECT_SOURCE_DIR}/testing/analysis_tests
                ${PROJECT_BINARY_DIR}/analysis_tests
        )
    endif(WIN32)
endif()

configure_file(
    ${PROJECT_SOURCE_DIR}//external_libs/SandiaDecay/sandia.decay.nocoinc.min.xml
    ${PROJECT_BINARY_DIR}/data/sandia.decay.xml
    COPYONLY
)
