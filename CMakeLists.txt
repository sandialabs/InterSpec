cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set( CMAKE_POSITION_INDEPENDENT_CODE ON CACHE BOOL "Set PIC to make InterSpec static lib" )


if( InterSpec_FETCH_DEPENDENCIES )
  include( cmake/FetchInterSpecDeps.txt )
endif( InterSpec_FETCH_DEPENDENCIES )


project(InterSpec VERSION 1.10.0)

# If we are the top-level project, we'll build the local-server exe
get_directory_property(hasParent PARENT_DIRECTORY)


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

if( hasParent )
  option( BUILD_AS_LOCAL_SERVER "Build for use as a server for localhost" OFF )
else( hasParent )
  option( BUILD_AS_LOCAL_SERVER "Build for use as a server for localhost" ON )
endif( hasParent )


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
    USE_SPECRUM_FILE_QUERY_WIDGET
    "Enables compilation of the spectrum file query tool (not for web deployment)"
    OFF
)
option(
    USE_SEARCH_MODE_3D_CHART
    "Enable 3D chart or searchmode/passthrough data.  Requires Wt >=3.3.4 compiled with WebGL/OpenGL support"
    ON
)
option(USE_TERMINAL_WIDGET "Enables use of terminal widget" ON)

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

if(MSVC)
  option(${PROJECT_NAME}_USE_MSVC_MultiThreadDLL "Use dynamically-link runtime library." OFF)

  if( ${PROJECT_NAME}_USE_MSVC_MultiThreadDLL)
    set(Boost_USE_STATIC_RUNTIME OFF)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
  else()
    set(Boost_USE_STATIC_RUNTIME ON)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  endif()
  message(STATUS "Using CMAKE_MSVC_RUNTIME_LIBRARY ${CMAKE_MSVC_RUNTIME_LIBRARY}")
endif()


set(MAX_SPECTRUM_MEMMORY_SIZE_MB 256
    CACHE STRING
    "Amount of memory to allow spectra to take up before trying to offload them onto disk when not in use"
)

set(GOOGLE_MAPS_KEY "" CACHE STRING "Google maps api key.")

set(GUI_TYPE "")
set( INTERSPEC_LIB_TYPE "STATIC" )
#set( INTERSPEC_LIB_TYPE "SHARED" )


if(WIN32)
#    set(EXECUTABLE_OUTPUT_PATH . CAxCHE PATH "Path to executables" FORCE)
    add_definitions(-D _WIN32_WINNT=0x0601 -D WINVER=0x0601 -D _SCL_SECURE_NO_WARNINGS ) #0x0601==Win7, 0x0501=WinXP (Wt CMakeLists.txt uses 0x0501x)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP /wd4996 /wd4267" )
else(WIN32)
#    set(EXECUTABLE_OUTPUT_PATH bin CACHE PATH "Path to executables" FORCE)
endif(WIN32)


SET(MYSQL_DATABASE_TO_USE "dev" CACHE STRING "Which mysql database to use (if applicable): dev/qc/prod")


# I think requiring c++14 will require GCC>=5, but we'll double check anyway
if( CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0 )
  message( fatal, "GCC less than 5 not supported." )
endif()


if(WIN32 OR BUILD_AS_ELECTRON_APP)
    set(TRY_TO_STATIC_LINK "ON")
endif(WIN32 OR BUILD_AS_ELECTRON_APP)

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


if(BUILD_AS_ELECTRON_APP)
    set(BUILD_AS_LOCAL_SERVER OFF)
    set(ALLOW_URL_TO_FILESYSTEM_MAP ON)
    set(USE_OSX_NATIVE_MENU OFF)
    set(USE_SPECRUM_FILE_QUERY_WIDGET ON)
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
endif(BUILD_FOR_WEB_DEPLOYMENT)

if(BUILD_AS_LOCAL_SERVER)
    set(PERFORM_DEVELOPER_CHECKS ON)
endif(BUILD_AS_LOCAL_SERVER)

#check to make sure exactly one target is specified
set(NUM_TARGETS_SPECIFIED 0)

if( BUILD_AS_OSX_APP )
  MATH( EXPR NUM_TARGETS_SPECIFIED "${NUM_TARGETS_SPECIFIED}+1" )
endif()

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
    if(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif(WIN32)
else(TRY_TO_STATIC_LINK)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif(TRY_TO_STATIC_LINK)

INCLUDE(cmake/FindWt.cmake)

find_package( ZLIB REQUIRED )
find_package( Threads REQUIRED )

#set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

find_package(Boost 1.68 REQUIRED COMPONENTS thread date_time system filesystem program_options regex random )



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
    src/ShieldingSelect.cpp
    src/GammaXsGui.cpp
    src/ReferenceLineInfo.cpp
    src/ReferencePhotopeakDisplay.cpp
    src/IsotopeNameFilterModel.cpp
    src/ReactionGamma.cpp
    src/IsotopeSearchByEnergy.cpp
    src/IsotopeSearchByEnergyModel.cpp
    src/SpectrumChart.cpp
    src/SpectrumDataModel.cpp
    src/HelpSystem.cpp
    src/InterSpec.cpp
    src/PopupDiv.cpp
    src/SpecMeas.cpp
    src/PeakFit.cpp
    src/PeakFitUtils.cpp
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
    src/DrfSelect.cpp
    src/MakeDrf.cpp
    src/MakeDrfSrcDef.cpp
    src/MakeDrfChart.cpp
    src/MakeDrfFit.cpp
    src/FluxTool.cpp
    src/NativeFloatSpinBox.cpp
    src/SwitchCheckbox.cpp
    src/EnergyCalTool.cpp
    src/EnergyCal.cpp
    src/EnergyCalGraphical.cpp
    src/EnergyCalAddActions.cpp
    src/EnergyCalPreserveWindow.cpp
    src/EnergyCalMultiFile.cpp
    src/SimpleDialog.cpp
    src/FeatureMarkerWidget.cpp
    src/ShowRiidInstrumentsAna.cpp
    src/D3SpectrumDisplayDiv.cpp
    external_libs/SpecUtils/d3_resources/SpectrumChartD3.js
    external_libs/SpecUtils/d3_resources/SpectrumChartD3.css
    external_libs/SpecUtils/d3_resources/SpectrumChartD3StandAlone.css
    src/D3TimeChart.cpp
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
    InterSpec/ShieldingSelect.h
    InterSpec/GammaXsGui.h
    InterSpec/ReferenceLineInfo.h
    InterSpec/ReferencePhotopeakDisplay.h
    InterSpec/IsotopeNameFilterModel.h
    InterSpec/ReactionGamma.h
    InterSpec/IsotopeSearchByEnergy.h
    InterSpec/IsotopeSearchByEnergyModel.h
    InterSpec/SpectrumChart.h
    InterSpec/SpectrumDataModel.h
    InterSpec/HelpSystem.h
    InterSpec/InterSpec.h
    InterSpec/PopupDiv.h
    InterSpec/SpecMeas.h
    InterSpec/PeakFit.h
    InterSpec/PeakFitUtils.h
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
    InterSpec/DrfSelect.h
    InterSpec/MakeDrf.h
    InterSpec/MakeDrfSrcDef.h
    InterSpec/MakeDrfChart.h
    InterSpec/MakeDrfFit.h
    InterSpec/FluxTool.h
    InterSpec/NativeFloatSpinBox.h
    InterSpec/SwitchCheckbox.h
    InterSpec/EnergyCalTool.h
    InterSpec/EnergyCal.h
    InterSpec/EnergyCalGraphical.h
    InterSpec/EnergyCalAddActions.h
    InterSpec/EnergyCalPreserveWindow.h
    InterSpec/EnergyCalMultiFile.h
    InterSpec/SimpleDialog.h
    InterSpec/FeatureMarkerWidget.h
    InterSpec/ShowRiidInstrumentsAna.h
    InterSpec/D3SpectrumDisplayDiv.h
    InterSpec/D3TimeChart.h
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
set(SpecUtils_ENABLE_D3_CHART ON CACHE BOOL "")
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
elseif(ANDROID)
    # TODO: Need to move making Android shared library to its CMakeLists.txt
    #list(APPEND headers target/android/AndroidUtils.hpp)
    #add_library(
    #    InterSpecLib
    #    SHARED
    #    main.cpp
    #    ${sources}
    #    ${headers}
    #)
elseif(BUILD_AS_UNIT_TEST_SUITE)
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
else()
endif(IOS)


add_library(InterSpecLib ${INTERSPEC_LIB_TYPE} ${sources} ${headers})


# If we are the top-level project, we'll build the local-server exe
if( NOT hasParent)
  add_executable(InterSpecExe ${GUI_TYPE} main.cpp)
  target_link_libraries(InterSpecExe PUBLIC InterSpecLib)
  set_target_properties(InterSpecExe PROPERTIES OUTPUT_NAME "InterSpec")
endif()


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
    Boost::thread
    Boost::date_time
    Boost::system
    Boost::filesystem
    Boost::program_options
    Boost::regex
    Boost::random
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
elseif(BUILD_AS_ELECTRON_APP)
elseif(YES)
endif(IOS)

if(NOT IOS)
    target_link_libraries(
        InterSpecLib
        PUBLIC
        ${Wt_LIBRARIES}
        ${Wt_LIBRARY}
        ${Wt_HTTP_LIBRARY}
        ${Wt_DBO_LIBRARY}
        ${Wt_DBOSQLITE3_LIBRARY}
        
        Threads::Threads
        Boost::system
        Boost::regex
        Boost::thread
        Boost::filesystem
        Boost::program_options
        Boost::date_time
        ZLIB::ZLIB
    )
endif(NOT IOS)

if(INCLUDE_ANALYSIS_TEST_SUITE)
  target_link_libraries(InterSpecLib PUBLIC ${Wt_TEST_LIBRARY})
endif(INCLUDE_ANALYSIS_TEST_SUITE)

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


#if(SUPPORT_ZIPPED_SPECTRUM_FILES)
#    #CMake seems to fail to find zlib on ios/android (even using find_package(ZLIB REQUIRED)), but we know it will be in
#    #  the systems path anyway, so just blindly link to zlib
#    if(NOT (IOS OR ANDROID))
#        if(NOT ZLIB_FOUND)
#            MESSAGE(FATAL_ERROR       "ZLIB required to enable SUPPORT_ZIPPED_SPECTRUM_FILES, but it wasnt found")
#        else(NOT ZLIB_FOUND)
#            target_include_directories(InterSpecLib PUBLIC ${ZLIB_INCLUDE_DIRS})
#            target_link_libraries(InterSpecLib PUBLIC ${ZLIB_LIBRARIES})
#        endif(NOT ZLIB_FOUND)
#    else(NOT (IOS OR ANDROID))
#        target_link_libraries(InterSpecLib PUBLIC z)
#    endif(NOT (IOS OR ANDROID))
#endif(SUPPORT_ZIPPED_SPECTRUM_FILES)

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

if( WIN32 )
  # The following are only needed for Wt 3.7.1, not 3.3.4
  target_link_libraries( InterSpecLib PRIVATE opengl32.lib d2d1.lib windowscodecs.lib Dwrite.lib )
endif( WIN32 )



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
