 cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(InterSpec VERSION 1.0.14)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set( CMAKE_POSITION_INDEPENDENT_CODE ON CACHE BOOL "Set PIC to make InterSpec static lib" )

if( CONTAINER_LIBC OR CONTAINER_MUSL)
  set( BUILD_AS_LOCAL_SERVER OFF CACHE INTERNAL "" )
  set( CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
  set( USE_REL_ACT_TOOL ON CACHE INTERNAL "")
  set( BUILD_FOR_WEB_DEPLOYMENT ON CACHE INTERNAL "" )
  set( CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  #set( Boost_INCLUDE_DIR "build/_deps/boost-src/libs")
  set( ZLIB_INCLUDE_DIR:STRING "build/_deps/zlib-build")
  set( SpecUtils_PYTHON_BINDINGS ON )
  set( InterSpec_FETCH_DEPENDENCIES ON CACHE INTERNAL "" )
  if( CONTAINER_MUSL )
    # In Alpine, the hardening and link-time-optimization clash for `vsnprintf`, so we need to disable one
    # of these - and for serving to web, lto seems like the better one to disable.
    #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -U_FORTIFY_SOURCE" ) # or  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORTIFY_SOURCE=0")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
    set( CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++" CACHE INTERNAL "")
    set( CMAKE_FIND_LIBRARY_SUFFIXES ".a" CACHE INTERNAL "")
  endif()
endif()

if( InterSpec_FETCH_DEPENDENCIES )
  include( cmake/FetchInterSpecDeps.cmake )
endif( InterSpec_FETCH_DEPENDENCIES )

# If we are the top-level project, we'll build the local-server exe
get_directory_property(hasParent PARENT_DIRECTORY)

option(
    BUILD_AS_ELECTRON_APP
    "Configures the executable to be ran as a Electron app"
    OFF
)
option(
    BUILD_AS_WX_WIDGETS_APP
    "Configures executable to be ran as a wxWidgets-based app"
    OFF
)
option(
    BUILD_FOR_WEB_DEPLOYMENT
    "Selects options appropriate for the web in the future"
    OFF
)
option(
    BUILD_FOR_CONTAINER_LIBC
    "Configures the executable to be run as a Docker app in a libc-based container"
)
option(
    BUILD_FOR_CONTAINER_MUSL
    "Configures the executable to be run as a Docker app in a musl-based container"
)
if( hasParent )
  option( BUILD_AS_LOCAL_SERVER "Build for use as a server for localhost" OFF )
else( hasParent )
  option( BUILD_AS_LOCAL_SERVER "Build for use as a server for localhost" ON )
endif( hasParent )


option(
    INCLUDE_ANALYSIS_TEST_SUITE
    "Allow whether user can save and load test spectra"
    OFF
)
option(
    PERFORM_DEVELOPER_CHECKS
    "Performs additional computationally expensive tests during execution"
    OFF
)
option(
    USE_OSX_NATIVE_MENU
    "Mirrors WMenus with OSX native NSMenu (not fully implemented)"
    OFF
)


if( NOT DEFINED TRY_TO_STATIC_LINK )
  if(WIN32)
    set( TRY_TO_STATIC_LINK ON )
  else(WIN32)
    option(TRY_TO_STATIC_LINK "Try to link to static libs whenever possible" OFF)
  endif(WIN32)
endif( NOT DEFINED TRY_TO_STATIC_LINK )

option(
    DRAW_GAMMA_LINES_LOG_AND_LIN
    "Have gamma lines indicators be drawn in log mode when y-axis is log"
    OFF
)
option(USE_DB_TO_STORE_SPECTRA "Use the database to store spectra" ON)
option(
    USE_SPECRUM_FILE_QUERY_WIDGET
    "Enables compilation of the spectrum file query tool (not for web deployment)"
    ON
)
option(
    USE_SEARCH_MODE_3D_CHART
    "Enable 3D chart or searchmode/passthrough data.  Requires Wt >=3.3.4 compiled with WebGL/OpenGL support"
    ON
)
option(USE_TERMINAL_WIDGET "Enables use of terminal widget" ON)

option(USE_DETECTION_LIMIT_TOOL "Enables under-development detection limit tool" ON)

option(USE_QR_CODES "Enables generating QR codes for like DRFs" ON)

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
option(USE_GOOGLE_MAP "Use google maps widget" OFF)
set(GOOGLE_MAPS_KEY "" CACHE STRING "Google maps api key.")

option(USE_LEAFLET_MAP "Use leaflet-based maps widget" OFF)
set(LEAFLET_MAPS_KEY "n/a" CACHE STRING "Leaflet maps api key.")

option(USE_REMOTE_RID "Enables using remote RID tool" ON)
option(USE_REL_ACT_TOOL "Enables Relative Activity tool" ON)
option(USE_LLM_INTERFACE "Enables the AI/LLM interface and MCP server code" ON)



if(BUILD_FOR_WEB_DEPLOYMENT)
    set(TRY_TO_STATIC_LINK OFF)
    set(CMAKE_BUILD_TYPE:STRING "Release")
    set(USE_DB_TO_STORE_SPECTRA ON)
    set(PERFORM_DEVELOPER_CHECKS OFF CACHE INTERNAL "")
    set(SpecUtils_ENABLE_EQUALITY_CHECKS OFF CACHE INTERNAL "")
    set(USE_SPECRUM_FILE_QUERY_WIDGET OFF)
endif(BUILD_FOR_WEB_DEPLOYMENT)


if( IOS OR ANDROID OR BUILD_FOR_WEB_DEPLOYMENT )
  set(USE_BATCH_CLI_TOOLS OFF CACHE INTERNAL "")
  set(USE_BATCH_GUI_TOOLS OFF CACHE INTERNAL "")
else()
  option( USE_BATCH_CLI_TOOLS "Enables command line batch processing." ON )
  option( USE_BATCH_GUI_TOOLS "Enables GUI batch processing." ON )
endif()

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

# If we are creating the normal InterSpec executable, and the InterSpec_batch executable,
#  we will compile LibInterSpec as a shared library, so both exes can use it; otherwise
#  we will create a static library, and statically link it into our final executable.
#  (Electron doesnt need a separate exe for batch, and we dont care abour disk space for unit tests)
if( USE_BATCH_CLI_TOOLS AND NOT BUILD_AS_ELECTRON_APP AND NOT BUILD_AS_UNIT_TEST_SUITE )
  set( INTERSPEC_LIB_TYPE SHARED )
  set( BUILD_LibInterSpec_SHARED ON )
else()
  set( INTERSPEC_LIB_TYPE STATIC )
  set( BUILD_LibInterSpec_SHARED OFF )
endif()


SET(MYSQL_DATABASE_TO_USE "dev" CACHE STRING "Which mysql database to use (if applicable): dev/qc/prod")

# Make a cached variable available to parent cmake projects
set( INTERSPEC_VERSION "${PROJECT_VERSION}" CACHE STRING "InterSpec version" )

if(NOT CMAKE_BUILD_TYPE)
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
endif(NOT CMAKE_BUILD_TYPE)


if(BUILD_AS_LOCAL_SERVER)
    set(PERFORM_DEVELOPER_CHECKS ON CACHE INTERNAL "")
    set(SpecUtils_ENABLE_EQUALITY_CHECKS ON CACHE INTERNAL "")
endif(BUILD_AS_LOCAL_SERVER)


if(BUILD_FOR_WEB_DEPLOYMENT AND BUILD_AS_LOCAL_SERVER )
  message( FATAL_ERROR "You can not build for both local webserver and web deployment" )
endif()



#try to link to static libraries whenever possible
if( TRY_TO_STATIC_LINK )
    if(WIN32)
        SET( ZLIB_USE_STATIC_LIBS "ON" )
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
    else(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
    endif(WIN32)
else(TRY_TO_STATIC_LINK)
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif(TRY_TO_STATIC_LINK)

include( cmake/FindWt.cmake )

find_package( ZLIB REQUIRED )
find_package( Threads REQUIRED )

if( InterSpec_FETCH_DEPENDENCIES )
  set( Wt_INCLUDE_DIR "${wt_SOURCE_DIR}/src" )
else( InterSpec_FETCH_DEPENDENCIES )
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_MULTITHREADED ON)
  find_package(Boost 1.68 REQUIRED COMPONENTS thread date_time system filesystem program_options regex random )

  find_package( Eigen3 3.4 REQUIRED NO_MODULE )
  find_package( Ceres REQUIRED )
endif( InterSpec_FETCH_DEPENDENCIES )



set(sources
    src/InterSpecApp.cpp
    src/DecayDataBaseServer.cpp
    src/IsotopeSelectionAids.cpp
    src/IsotopeId.cpp
    src/MaterialDB.cpp
    src/PhysicalUnits.cpp
    src/PhysicalUnitsLocalized.cpp
    src/MassAttenuationTool.cpp
    src/DetectorPeakResponse.cpp
    src/GammaInteractionCalc.cpp
    src/ShieldingSourceDisplay.cpp
    src/ShieldingSourceDiagram.cpp
    src/ShieldingSourceFitCalc.cpp
    src/ShieldingSelect.cpp
    src/GammaXsGui.cpp
    src/ReferenceLineInfo.cpp
    src/ReferenceLinePredef.cpp
    src/ReferencePhotopeakDisplay.cpp
    src/IsotopeNameFilterModel.cpp
    src/ReactionGamma.cpp
    src/IsotopeSearchByEnergy.cpp
    src/IsotopeSearchByEnergyModel.cpp
    src/SpectrumChart.cpp
    src/SpectrumDataModel.cpp
    src/HelpSystem.cpp
    src/InjaLogDialog.cpp
    src/InterSpec.cpp
    src/InterSpecServer.cpp
    src/PopupDiv.cpp
    src/SpecMeas.cpp
    src/PeakFit.cpp
    src/PeakFitUtils.cpp
    src/PeakDef.cpp
    src/PeakDists.cpp
    src/XRayWidthServer.cpp
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
    src/FileDragUploadResource.cpp
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
    src/DrfChart.cpp
    src/MakeDrf.cpp
    src/MakeDrfSrcDef.cpp
    src/MakeDrfChart.cpp
    src/MakeDrfFit.cpp
    src/MakeFwhmForDrf.cpp
    src/FluxTool.cpp
    src/NativeFloatSpinBox.cpp
    src/SwitchCheckbox.cpp
    src/EnergyCalTool.cpp
    src/EnergyCal.cpp
    src/EnergyCalGraphical.cpp
    src/EnergyCalAddActions.cpp
    src/EnergyCalPreserveWindow.cpp
    src/EnergyCalMultiFile.cpp
    src/FitPeaksForNuclides.cpp
    src/FitPeaksForNuclidesGui.cpp
    src/SimpleDialog.cpp
    src/FeatureMarkerWidget.cpp
    src/ShowRiidInstrumentsAna.cpp
    src/D3SpectrumDisplayDiv.cpp
    external_libs/SpecUtils/d3_resources/SpectrumChartD3.js
    external_libs/SpecUtils/d3_resources/SpectrumChartD3.css
    src/D3TimeChart.cpp
    src/ZipArchive.cpp
    src/MultimediaDisplay.cpp
    src/MoreNuclideInfo.cpp
    src/MoreNuclideInfoDisplay.cpp
    src/UndoRedoManager.cpp
    src/AppUtils.cpp
    src/ExportSpecFile.cpp
    src/EnterAppUrlWindow.cpp
    src/ExternalRidResult.cpp
    src/NuclideSourceEnter.cpp
    src/AddNewPeakDialog.cpp
    src/UserPreferences.cpp
    src/RefSpectraWidget.cpp
    src/RefSpectraModel.cpp
    src/SimpleActivityCalc.cpp
    src/RefLineDynamic.cpp
    src/AnalystChecks.cpp
    src/BatchPeak.cpp
    src/BatchActivity.cpp
    src/BatchInfoLog.cpp
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
    InterSpec/PhysicalUnitsLocalized.h
    InterSpec/MassAttenuationTool.h
    InterSpec/DetectorPeakResponse.h
    InterSpec/GammaInteractionCalc.h
    InterSpec/ShieldingSourceDisplay.h
    InterSpec/ShieldingSourceDiagram.h
    InterSpec/ShieldingSourceFitCalc.h
    InterSpec/ShieldingSelect.h
    InterSpec/GammaXsGui.h
    InterSpec/ReferenceLineInfo.h
    InterSpec/ReferenceLinePredef.h
    InterSpec/ReferencePhotopeakDisplay.h
    InterSpec/IsotopeNameFilterModel.h
    InterSpec/ReactionGamma.h
    InterSpec/IsotopeSearchByEnergy.h
    InterSpec/IsotopeSearchByEnergyModel.h
    InterSpec/SpectrumChart.h
    InterSpec/SpectrumDataModel.h
    InterSpec/HelpSystem.h
    InterSpec/InjaLogDialog.h
    InterSpec/InterSpec.h
    InterSpec/InterSpecServer.h
    InterSpec/PopupDiv.h
    InterSpec/SpecMeas.h
    InterSpec/PeakFit.h
    InterSpec/PeakFit_imp.hpp
    InterSpec/PeakFitUtils.h
    InterSpec/PeakDef.h
    InterSpec/PeakDists.h
    InterSpec/PeakDists_imp.hpp
    InterSpec/XRayWidthServer.h
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
    InterSpec/FileDragUploadResource.h
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
    InterSpec/DrfChart.h
    InterSpec/MakeDrf.h
    InterSpec/MakeDrfSrcDef.h
    InterSpec/MakeDrfChart.h
    InterSpec/MakeDrfFit.h
    InterSpec/MakeFwhmForDrf.h
    InterSpec/FluxTool.h
    InterSpec/NativeFloatSpinBox.h
    InterSpec/SwitchCheckbox.h
    InterSpec/EnergyCalTool.h
    InterSpec/EnergyCal.h
    InterSpec/EnergyCalGraphical.h
    InterSpec/EnergyCalAddActions.h
    InterSpec/EnergyCalPreserveWindow.h
    InterSpec/EnergyCalMultiFile.h
    InterSpec/FitPeaksForNuclides.h
    InterSpec/FitPeaksForNuclidesGui.h
    InterSpec/SimpleDialog.h
    InterSpec/FeatureMarkerWidget.h
    InterSpec/ShowRiidInstrumentsAna.h
    InterSpec/D3SpectrumDisplayDiv.h
    InterSpec/D3TimeChart.h
    InterSpec/ZipArchive.h
    InterSpec/MultimediaDisplay.h
    InterSpec/MoreNuclideInfo.h
    InterSpec/MoreNuclideInfoDisplay.h
    InterSpec/UndoRedoManager.h
    InterSpec/AppUtils.h
    InterSpec/BersteinPolynomial.hpp
    InterSpec/ExportSpecFile.h
    InterSpec/ExternalRidResult.h
    InterSpec/EnterAppUrlWindow.h
    InterSpec/NuclideSourceEnter.h
    InterSpec/AddNewPeakDialog.h
    InterSpec/UserPreferences.h
    InterSpec/XmlUtils.hpp
    InterSpec/RefSpectraWidget.h
    InterSpec/RefSpectraModel.h
    InterSpec/SimpleActivityCalc.h
    InterSpec/RefLineDynamic.h
    InterSpec/AnalystChecks.h
    InterSpec/BatchPeak.h
    InterSpec/BatchActivity.h
    InterSpec/BatchInfoLog.h
)

if(USE_DB_TO_STORE_SPECTRA)
    list(APPEND sources src/DbFileBrowser.cpp)
    list(APPEND headers InterSpec/DbFileBrowser.h)
endif(USE_DB_TO_STORE_SPECTRA)


if(INCLUDE_ANALYSIS_TEST_SUITE)
    set(TEST_SUITE_BASE_DIR "." CACHE STRING "Directory containing a 'analysis_tests' directory containing N42 analyst test files.")
    list(APPEND sources src/SpectrumViewerTester.cpp)
    list(APPEND headers InterSpec/SpectrumViewerTester.h)
endif(INCLUDE_ANALYSIS_TEST_SUITE)


#if(IOS)
#    list(APPEND sources target/ios/InterSpec/FileHandling.mm)
#    list(APPEND headers target/ios/InterSpec/FileHandling.h)
#endif(IOS)

if(USE_GOOGLE_MAP)
    list(APPEND sources src/GoogleMap.cpp)
    list(APPEND headers InterSpec/GoogleMap.h)
endif(USE_GOOGLE_MAP)

if(USE_LEAFLET_MAP)
    list(APPEND sources src/LeafletRadMap.cpp)
    list(APPEND headers InterSpec/LeafletRadMap.h)
    list(APPEND OTHER_InterSpec_SUPPORT_FILES InterSpec_resources/LeafletRadMap.js)
endif(USE_LEAFLET_MAP)

if(USE_SEARCH_MODE_3D_CHART)
    list(APPEND sources src/SearchMode3DChart.cpp src/SearchMode3DDataModel.cpp)
    list( APPEND headers InterSpec/SearchMode3DChart.h InterSpec/SearchMode3DDataModel.h )
endif(USE_SEARCH_MODE_3D_CHART)

if(USE_TERMINAL_WIDGET)
    list( APPEND sources src/TerminalModel.cpp src/TerminalWidget.cpp js/TerminalWidget.js )
    list( APPEND headers InterSpec/TerminalModel.h InterSpec/TerminalWidget.h )
endif(USE_TERMINAL_WIDGET)

if(USE_DETECTION_LIMIT_TOOL)
    list( APPEND sources src/DetectionLimitCalc.cpp src/DetectionLimitTool.cpp src/DetectionLimitSimple.cpp )
    list( APPEND headers InterSpec/DetectionLimitCalc.h InterSpec/DetectionLimitTool.h InterSpec/DetectionLimitSimple.h )
endif(USE_DETECTION_LIMIT_TOOL)

if(USE_QR_CODES)
  list( APPEND sources external_libs/QR-Code-generator/cpp/qrcodegen.cpp 
        src/QrCode.cpp 
        src/QRSpectrum.cpp
  )

  list( APPEND headers external_libs/QR-Code-generator/cpp/qrcodegen.hpp 
        InterSpec/QrCode.h 
        InterSpec/QRSpectrum.h 
  )
endif(USE_QR_CODES)

if(USE_SPECRUM_FILE_QUERY_WIDGET)
  if( BUILD_FOR_WEB_DEPLOYMENT )
    message( FATAL_ERROR "The spectrum file query widget can not be enabled for web deployments." )
  endif( BUILD_FOR_WEB_DEPLOYMENT )
        
  list( APPEND sources src/SpecFileQuery.cpp src/SpecFileQueryWidget.cpp src/SpecFileQueryDbCache.cpp )
  list( APPEND headers InterSpec/SpecFileQuery.h InterSpec/SpecFileQueryWidget.h InterSpec/SpecFileQueryDbCache.h )
  list( APPEND sources js/SpecFileQueryWidget.js )
endif(USE_SPECRUM_FILE_QUERY_WIDGET)


if( USE_REMOTE_RID )
  list( APPEND sources src/RemoteRid.cpp  )
  list( APPEND headers InterSpec/RemoteRid.h )
endif( USE_REMOTE_RID )


if( USE_REL_ACT_TOOL )
  list( APPEND sources
        src/PeakFitLM.cpp
        src/RelActCalc.cpp
        src/RelActCalcAuto.cpp
        src/RelActCalcManual.cpp
        src/RelActAutoGui.cpp
        src/RelActManualGui.cpp
        src/RelActTxtResults.cpp
        src/RelEffChart.cpp
        src/RelEffShieldWidget.cpp
        src/RelActAutoGuiRelEffOptions.cpp
        src/RelActAutoGuiEnergyRange.cpp
        src/RelActAutoGuiNuclide.cpp
        src/RelActAutoGuiFreePeak.cpp
        src/RelActAutoReport.cpp
        src/RelActAutoDev.cpp
        src/ShieldingSourceFitPlot.cpp
  )
  
  list( APPEND headers
        InterSpec/PeakFitLM.h
        InterSpec/RelActCalc.h
        InterSpec/RelActCalc_imp.hpp
        InterSpec/RelActCalcAuto.h
        InterSpec/RelActCalcAuto_imp.hpp
        InterSpec/RelActCalcAuto_EnergyCal_imp.hpp
        InterSpec/RelActCalcManual.h
        InterSpec/RelActAutoGui.h
        InterSpec/RelActManualGui.h
        InterSpec/RelActTxtResults.h
        InterSpec/RelEffChart.h
        InterSpec/RelEffShieldWidget.h
        InterSpec/ShieldingSourceFitPlot.h
        InterSpec/RelActAutoGuiRelEffOptions.h
        InterSpec/RelActAutoGuiEnergyRange.h
        InterSpec/RelActAutoGuiNuclide.h
        InterSpec/RelActAutoGuiFreePeak.h
        InterSpec/RelActAutoReport.h
        InterSpec/RelActAutoDev.h
  )
endif( USE_REL_ACT_TOOL )


if( USE_BATCH_CLI_TOOLS )
  list( APPEND sources src/BatchCommandLine.cpp )
  list( APPEND headers InterSpec/BatchCommandLine.h )
endif( USE_BATCH_CLI_TOOLS )


if( USE_BATCH_GUI_TOOLS )
  list( APPEND sources
        src/BatchGuiWidget.cpp
        src/BatchGuiAnaWidget.cpp
        src/BatchGuiInputFile.cpp
  )
  
  list( APPEND headers
        InterSpec/BatchGuiWidget.h
        InterSpec/BatchGuiAnaWidget.h
        InterSpec/BatchGuiInputFile.h
  )
endif( USE_BATCH_GUI_TOOLS )

if( (NOT IOS) AND (NOT ANDROID) AND (NOT BUILD_FOR_WEB_DEPLOYMENT) )
  list( APPEND sources src/DirectorySelector.cpp )
  list( APPEND headers InterSpec/DirectorySelector.h )
endif( (NOT IOS) AND (NOT ANDROID) AND (NOT BUILD_FOR_WEB_DEPLOYMENT) )


if( BUILD_AS_OSX_APP )
  list( APPEND sources target/osx/macOsUtils.mm target/osx/NativeMenu.mm )
  list( APPEND headers target/osx/macOsUtils.h target/osx/NativeMenu.h )
endif( BUILD_AS_OSX_APP )

if( USE_LLM_INTERFACE )
  list( APPEND sources src/LlmMcpResource.cpp src/LlmInterface.cpp src/LlmToolRegistry.cpp src/LlmIsotopicsTool.cpp src/LlmActivityFitTool.cpp src/LlmRelActManualTool.cpp src/LlmConfig.cpp src/LlmConversationHistory.cpp src/LlmToolGui.cpp src/LlmInteractionDisplay.cpp )
  list( APPEND headers InterSpec/LlmMcpResource.h InterSpec/LlmInterface.h InterSpec/LlmToolRegistry.h InterSpec/LlmIsotopicsTool.h InterSpec/LlmActivityFitTool.h InterSpec/LlmRelActManualTool.h InterSpec/LlmConfig.h InterSpec/LlmConversationHistory.h InterSpec/LlmToolGui.h InterSpec/LlmInteractionDisplay.h )
  list( APPEND OTHER_InterSpec_SUPPORT_FILES data/llm_config.xml data/llm_tools_config.xml data/llm_agents.xml )
  
  if( BUILD_FOR_WEB_DEPLOYMENT OR IOS OR ANDROID )
    MESSAGE(FATAL_ERROR "MCP server can not be enabled for web or mobile builds")
  endif( BUILD_FOR_WEB_DEPLOYMENT OR IOS OR ANDROID )
endif( USE_LLM_INTERFACE )

#Copy D3 resources into InterSpec_resources directory at compile time
include(cmake/DeployJsAndCss.cmake)

list(APPEND sources external_libs/SpecUtils/d3_resources/d3.v3.min.js)
deploy_js_resource("${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/d3_resources/d3.v3.min.js"
                    "${CMAKE_CURRENT_BINARY_DIR}/InterSpec_resources/d3.v3.min.js"
)



list( APPEND OTHER_InterSpec_SUPPORT_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/D3TimeChart.js
  ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/ShieldingSourceFitPlot.js
  ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/ShieldingSourceFitPlot.css
  ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/Shielding3DView.js
  ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/Shielding3DView.css
  ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/Shielding2DView.js
  ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/Shielding2DView.css
)

deploy_js_resource("${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/d3_resources/SpectrumChartD3.js"
         "${CMAKE_CURRENT_BINARY_DIR}/InterSpec_resources/SpectrumChartD3.js"
         
)

set( CSS_RESOURCES_TO_COPY "SpectrumChartD3.css" )
    
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
set(SpecUtils_NO_BOOST_LIB OFF
    CACHE BOOL "Use boost for threading primitives in SpecUtils"
)

if( APPLE OR ANDROID )
  set(SpecUtils_USE_WT_THREADPOOL OFF CACHE BOOL "")
else( APPLE OR ANDROID  )
  set(SpecUtils_USE_WT_THREADPOOL ON CACHE BOOL "")
endif( APPLE OR ANDROID  )
set(SpecUtils_EXTERNALLY_DEFINED_LOG_MESSAGE ON CACHE BOOL "")
set(SpecUtils_BUILD_UNIT_TESTS OFF CACHE BOOL "")
set(SpecUtils_BUILD_REGRESSION_TEST OFF CACHE BOOL "")
set(SpecUtils_PYTHON_BINDINGS OFF CACHE BOOL "")
set(SpecUtils_JAVA_SWIG OFF CACHE BOOL "")
set(SpecUtils_ENABLE_D3_CHART ON CACHE BOOL "")
set(SpecUtils_D3_SUPPORT_FILE_STATIC OFF CACHE BOOL "")
set(SpecUtils_D3_SCRIPTS_RUNTIME_DIR "InterSpec_resources" CACHE STRING "")
set(SpecUtils_ENABLE_URI_SPECTRA ${USE_QR_CODES} CACHE BOOL "")
set(SpecUtils_FLT_PARSE_METHOD "boost" CACHE STRING "")

#add in all the subdirectories the project will always use
add_subdirectory(external_libs/Minuit2)
add_subdirectory(external_libs/Cuba-3.0)
add_subdirectory(external_libs/SpecUtils)
add_subdirectory(external_libs/SandiaDecay)
add_subdirectory(external_libs/muparserx-4.0.7)
add_subdirectory(external_libs/VoigtDistribution)

add_library(InterSpecLib ${INTERSPEC_LIB_TYPE} ${sources} ${headers} )

# If we are the top-level project, we'll build the local-server exe
if( NOT hasParent)
  add_executable(InterSpecExe main.cpp)
  target_link_libraries(InterSpecExe PUBLIC InterSpecLib)
  set_target_properties(InterSpecExe PROPERTIES OUTPUT_NAME "InterSpec")
  # Tell the runtime linker to search ../lib beside the executable
  if(UNIX AND NOT APPLE)           # $ORIGIN is a Linux/BSD feature
    set_target_properties(InterSpecExe PROPERTIES
        BUILD_RPATH  "$ORIGIN/../lib"
        INSTALL_RPATH "$ORIGIN/../lib")
    # makes variables such as $ORIGIN legal in BUILD_RPATH
    set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
  endif()
endif()


set_target_properties(InterSpecLib PROPERTIES PREFIX "")
set_target_properties(InterSpecLib PROPERTIES OUTPUT_NAME "InterSpec")

if((UNIX AND NOT APPLE) AND (BUILD_FOR_WEB_DEPLOYMENT OR BUILD_AS_LOCAL_SERVER) )           # $ORIGIN is a Linux/BSD feature
  set_target_properties(InterSpecExe PROPERTIES
      BUILD_RPATH  "$ORIGIN/../lib"
      INSTALL_RPATH "$ORIGIN/../lib")
  # makes variables such as $ORIGIN legal in BUILD_RPATH
  set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
endif()

# Turn off some compiler warnings
target_compile_options(InterSpecLib PRIVATE
     $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
          -Wno-inconsistent-missing-override>  # Turn off a bunch of warnings like: WTreeView:399:20: warning: 'headerWidget' overrides a member function but is not marked 'override' [-Winconsistent-missing-override]
     $<$<CXX_COMPILER_ID:MSVC>:
          >)

# Turn off "note: The practice of declaring the Bind placeholders (_1, _2, ...) in the global namespace is deprecated...."
# warnings (these come from Wt I think, not our code)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -D BOOST_BIND_GLOBAL_PLACEHOLDERS")
          
if(WIN32)
  #0x0601==Win7, 0x0501=WinXP (Wt CMakeLists.txt uses 0x0501x)
  target_compile_definitions( InterSpecLib PUBLIC /bigobj -D _WIN32_WINNT=0x0601 -D WINVER=0x0601 -D _SCL_SECURE_NO_WARNINGS )
  target_compile_definitions( InterSpecLib PRIVATE InterSpec_EXPORTS )
  add_definitions( /bigobj )
        
  if( INTERSPEC_LIB_TYPE STREQUAL "SHARED" )
    # We get error compiling multi-processor debug build of shared LibInterSpec, presumably due to some 
    #  file locking issues, so only enable parallel builds for release, if we are building shared LibInterSpec
    target_compile_options( InterSpecLib PUBLIC $<$<CONFIG:Release>:/MP> )
  else()
    add_definitions( /MP )
  endif()
endif(WIN32)

target_link_libraries(
    InterSpecLib
    PUBLIC
    Threads::Threads
    SpecUtils
    SandiaDecay
    Minuit2::Minuit2
    Cuba-3.0
    muparserx
    Boost::thread
    Boost::date_time
    Boost::system
    Boost::filesystem
    Boost::program_options
    Boost::regex
    Boost::random
    ZLIB::ZLIB
    Eigen3::Eigen
    Ceres::ceres
    ${CMAKE_DL_LIBS}
    VoigtDistribution
)

if( USE_SPECRUM_FILE_QUERY_WIDGET )
    add_subdirectory(external_libs/pugixml-1.9)
    target_link_libraries(InterSpecLib PUBLIC pugixml)
endif( USE_SPECRUM_FILE_QUERY_WIDGET )


if(USE_MYSQL_DB EQUAL USE_SQLITE3_DB)
    MESSAGE(FATAL_ERROR "Exactly one type of database (MySQL or SQLITE3) must be chosen")
endif()

set( RAPID_XML_INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_libs/SpecUtils/3rdparty )

target_include_directories( InterSpecLib
  PUBLIC
    ${Wt_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
    ${RAPID_XML_INC_DIR}
    external_libs
)


if( InterSpec_FETCH_DEPENDENCIES )
  target_link_libraries( InterSpecLib
    PUBLIC
      Boost::math Boost::scope_exit Boost::uuid Boost::iostreams Boost::assign Boost::multi_array Boost::crc Boost::process Boost::crc
      wt
      wthttp
      wtdbosqlite3
      wtdbo
  )

  target_include_directories( InterSpecLib
    PUBLIC
      ${wt_BINARY_DIR}
      ${boost_SOURCE_DIR}/libs/numeric/ublas/include
  )
  
  # Need to fix copying the "resources" directory and configs and stuff probably
endif( InterSpec_FETCH_DEPENDENCIES )

# We will explicitly set the Wt include directory, marked as SYSTEM, to
#  avoid a ton of warning messages from those headers (if you take this
#  next line out, things will still work just fine, just with more warnings)
target_include_directories( InterSpecLib PUBLIC SYSTEM ${Wt_INCLUDE_DIR} )

if(WIN32)
    if( NOT InterSpec_FETCH_DEPENDENCIES )
      if( BUILD_AS_UNIT_TEST_SUITE )
        target_link_libraries(InterSpecLib PUBLIC  debug "${Wt_TEST_DEBUG_LIBRARY}" optimized "${Wt_TEST_LIBRARY}" )
      else( BUILD_AS_UNIT_TEST_SUITE )
        target_link_libraries(InterSpecLib PUBLIC  debug "${Wt_HTTP_DEBUG_LIBRARY}" optimized "${Wt_HTTP_LIBRARY}" )
      endif( BUILD_AS_UNIT_TEST_SUITE )
      
      target_link_libraries(InterSpecLib PUBLIC
        debug ${Wt_DEBUG_LIBRARY} ${Wt_DBO_DEBUG_LIBRARY} ${Wt_DBOSQLITE3_DEBUG_LIBRARY}
        optimized ${Wt_LIBRARY} ${Wt_DBO_LIBRARY} ${Wt_DBOSQLITE3_LIBRARY}
      )
    endif( NOT InterSpec_FETCH_DEPENDENCIES )
    
    target_link_directories(InterSpecLib PUBLIC ${Boost_LIBRARY_DIR})
else(WIN32)
  if( BUILD_AS_UNIT_TEST_SUITE )
    target_link_libraries(InterSpecLib PUBLIC ${Wt_TEST_LIBRARY})
  else( BUILD_AS_UNIT_TEST_SUITE )
    target_link_libraries(InterSpecLib PUBLIC ${Wt_HTTP_LIBRARY})
  endif( BUILD_AS_UNIT_TEST_SUITE )
  
  target_link_libraries(InterSpecLib PUBLIC ${Wt_LIBRARY} ${Wt_DBO_LIBRARY} ${Wt_DBOSQLITE3_LIBRARY} )
endif(WIN32)



target_link_libraries(
    InterSpecLib
    PUBLIC
    Threads::Threads
    SpecUtils
    SandiaDecay
    Minuit2::Minuit2
    Cuba-3.0
    muparserx
    Boost::thread
    Boost::date_time
    Boost::system
    Boost::filesystem
    Boost::program_options
    Boost::regex
    Boost::random
    ZLIB::ZLIB
    Eigen3::Eigen
    Ceres::ceres
    ${CMAKE_DL_LIBS}
)


if( BUILD_AS_OSX_APP )
  FIND_LIBRARY(COCOA_LIBRARY Cocoa)
  FIND_LIBRARY(COREDATA_LIBRARY CoreData)
  FIND_LIBRARY(APPKIT_LIBRARY AppKit)
  FIND_LIBRARY(FOUNDATION_LIBRARY Foundation)

  MARK_AS_ADVANCED( COCOA_LIBRARY COREDATA_LIBRARY APPKIT_LIBRARY FOUNDATION_LIBRARY )

  target_link_libraries( InterSpecLib PUBLIC
    ${COCOA_LIBRARY}
    ${APPKIT_LIBRARY}
    ${COREDATA_LIBRARY}
    ${FOUNDATION_LIBRARY}
  )
  
  if( INTERSPEC_LIB_TYPE STREQUAL "SHARED" )
    set_target_properties( InterSpecLib PROPERTIES
      FRAMEWORK TRUE
      FRAMEWORK_VERSION A
      MACOSX_FRAMEWORK_IDENTIFIER gov.sandia.InterSpecFramework
      MACOSX_FRAMEWORK_SHORT_VERSION_STRING 1.0
      MACOSX_FRAMEWORK_BUNDLE_VERSION 1.0.0
      MACOSX_FRAMEWORK_COPYRIGHT "© 2024 Sandia National Laboratories"
      # XCODE_ATTRIBUTE_INSTALL_PATH @executable_path/../Frameworks/
      XCODE_ATTRIBUTE_SKIP_INSTALL "Yes"
    )
    
    set_target_properties( InterSpecLib PROPERTIES
      XCODE_ATTRIBUTE_PRODUCT_BUNDLE_IDENTIFIER gov.sandia.InterSpecFramework
    )
    
    set_target_properties( InterSpecLib PROPERTIES
      PUBLIC_HEADER "InterSpec_config.h"
    )
  endif( INTERSPEC_LIB_TYPE STREQUAL "SHARED" )
endif( BUILD_AS_OSX_APP )


if( BUILD_AS_LOCAL_SERVER AND (CMAKE_GENERATOR STREQUAL "Xcode") )
    # Setting Xcode working directory requires cmake 3.17
    if( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17.0" )
        set_target_properties(InterSpecExe PROPERTIES
            XCODE_GENERATE_SCHEME TRUE
            XCODE_SCHEME_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        )
    endif( ${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.17.0" )
    
    set(CMAKE_XCODE_ATTRIBUTE_CODE_SIGNING_ALLOWED "NO") #macOS 26 requires this, or else we get some error for development
endif( BUILD_AS_LOCAL_SERVER AND (CMAKE_GENERATOR STREQUAL "Xcode") )

if(MSVC AND (NOT hasParent))
  set_target_properties( InterSpecExe PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif(MSVC AND (NOT hasParent))

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
        message(FATAL_ERROR "Could not find libpng - requred for libharu")
    else()
        target_link_libraries(InterSpecLib PUBLIC ${PNG_LIBRARY})
    endif()
endif()

#bellow is what I need to get things working 20131120 for statically linked
if(TRY_TO_STATIC_LINK)
    if(APPLE AND NOT IOS)
        find_library( CARBON_LIBRARY Carbon REQUIRED )
        find_library( FOUNDATION_LIBRARY Foundation REQUIRED )
        mark_as_advanced( CARBON_LIBRARY FOUNDATION_LIBRARY )
        target_link_libraries( InterSpecLib PUBLIC ${CARBON_LIBRARY} ${FOUNDATION_LIBRARY} )

        if(USE_OSX_NATIVE_MENU)
          find_library( APP_KIT_LIBRARY AppKit REQUIRED )
          mark_as_advanced( APP_KIT_LIBRARY )
          target_link_libraries( InterSpecLib PUBLIC ${APP_KIT_LIBRARY} )
        endif(USE_OSX_NATIVE_MENU)
    endif(APPLE AND NOT IOS)

    if(USE_MYSQL_DB)
        include(cmake/WtFindMysql.txt)
        target_link_libraries(InterSpecLib PUBLIC ${MYSQL_LIBRARIES})
    endif(USE_MYSQL_DB)

    if(Wt_FCGI_LIBRARY)
        include(cmake/WtFindFcgi.txt)
        target_link_libraries(InterSpecLib PUBLIC ${FCGI_LIBRARIES})
    endif(Wt_FCGI_LIBRARY)

    if(("${CMAKE_SYSTEM}" MATCHES "Linux") AND (NOT ANDROID))
        find_library( RT_LIB rt REQUIRES )
        target_link_libraries(InterSpecLib PUBLIC ${RT_LIB})

        #target_link_libraries( InterSpecExe -static-libgcc -static-libstdc++)
        #Next line statically links to pthread, copying whole archive into exe
        #set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-static -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive")
    endif(("${CMAKE_SYSTEM}" MATCHES "Linux") AND (NOT ANDROID))
endif(TRY_TO_STATIC_LINK)


if(Wt_DBOSQLITE3_LIBRARY)
    find_file( HAS_WT_DBO_SQLITE3 "Sqlite3" PATHS "${Wt_INCLUDE_DIR}/Wt/Dbo/backend/" NO_DEFAULT_PATH )
endif(Wt_DBOSQLITE3_LIBRARY)

if(Wt_MYSQL_LIBRARY)
    find_file( HAS_WT_DBO_MYSQL "MySQL" PATHS "${Wt_INCLUDE_DIR}/Wt/Dbo/backend/" NO_DEFAULT_PATH )
endif(Wt_MYSQL_LIBRARY)


if(Wt_MYSQL_LIBRARY AND BUILD_FOR_WEB_DEPLOYMENT AND USE_MYSQL_DB)
    set(DATABASE_PASSWORD_FILE "/path/to/passwords/databases.xml"
        CACHE STRING
        "Location of the XML passwords file for connecting to MySQL database"
    )
endif(Wt_MYSQL_LIBRARY AND BUILD_FOR_WEB_DEPLOYMENT AND USE_MYSQL_DB)

if(USE_MYSQL_DB AND NOT HAS_WT_DBO_MYSQL)
    MESSAGE(FATAL_ERROR   "You do not have the MySQL Wt Dbo backend installed, but you have asked to use it")
endif(USE_MYSQL_DB AND NOT HAS_WT_DBO_MYSQL)

if(USE_SQLITE3_DB AND NOT HAS_WT_DBO_SQLITE3 AND NOT WIN32 AND NOT InterSpec_FETCH_DEPENDENCIES )
    MESSAGE(FATAL_ERROR   "You do not have the SQLITE3 Wt Dbo backend installed, but you have asked to use it")
endif(USE_SQLITE3_DB AND NOT HAS_WT_DBO_SQLITE3 AND NOT WIN32 AND NOT InterSpec_FETCH_DEPENDENCIES )

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec/InterSpec_config.h.in
    ${CMAKE_BINARY_DIR}/InterSpec_config.h
)

#if(ANDROID)
#    target_include_directories(InterSpecLib PUBLIC target)
#endif(ANDROID)

if( WIN32 )
  # The following are only needed for Wt 3.7.1, not 3.3.4
  target_link_libraries( InterSpecLib PUBLIC opengl32.lib d2d1.lib windowscodecs.lib Dwrite.lib )
endif( WIN32 )

# after the target_link_libraries() for InterSpecLib / InterSpecExe
if (CMAKE_CXX_STANDARD LESS 17)
    # GCC < 9 still needs libstdc++fs for <experimental/filesystem>
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_link_libraries(InterSpecLib PUBLIC stdc++fs)
    endif()
endif()

#if( ${PERFORM_DEVELOPER_CHECKS} EQUAL ${SpecUtils_ENABLE_EQUALITY_CHECKS} )
#  MESSAGE(FATAL_ERROR "PERFORM_DEVELOPER_CHECKS (${PERFORM_DEVELOPER_CHECKS}) must be equal to SpecUtils_ENABLE_EQUALITY_CHECKS (${SpecUtils_ENABLE_EQUALITY_CHECKS})")
#endif()


# I cant seems to get more complex regexs to work, so we'll use many simple ones.
FILE(GLOB JS_BASE_FILES InterSpec_resources/*.js )
FILE(GLOB_RECURSE CSS_BASE_FILES LIST_DIRECTORIES TRUE InterSpec_resources/*.css )
list(FILTER CSS_BASE_FILES EXCLUDE REGEX ".*SpectrumChartD3.*css")
list(FILTER CSS_BASE_FILES EXCLUDE REGEX ".min.css")
FILE(GLOB STATIC_XML_FILES InterSpec_resources/static_text/*.xml )
FILE(GLOB STATIC_HTML_TMPLTS_FILES InterSpec_resources/static_text/*.html )
FILE(GLOB STATIC_ADD_XML_FILES InterSpec_resources/static_text/*/*.xml )
FILE(GLOB STATIC_ADD_CSS_FILES InterSpec_resources/static_text/*/*.css )

list( APPEND OTHER_InterSpec_SUPPORT_FILES
     ${JS_BASE_FILES}
     ${CSS_BASE_FILES}
     ${STATIC_XML_FILES}
     ${STATIC_ADD_XML_FILES}
     ${STATIC_ADD_CSS_FILES}
     ${STATIC_HTML_TMPLTS_FILES}
     ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec_resources/static_text/help.json
     ${CMAKE_CURRENT_SOURCE_DIR}/README.md
  )

add_custom_target(interspec_resources SOURCES ${OTHER_InterSpec_SUPPORT_FILES})


# Make a target for the XML files that hold the (internationalized) text for the app
FILE(GLOB APP_TEXT_FILES InterSpec_resources/app_text/*.xml )
add_custom_target( interspec_app_text SOURCES ${APP_TEXT_FILES} )


#Copy files to the install directory
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
    if(BUILD_FOR_WEB_DEPLOYMENT)
      file(COPY ${WT_RESOURCES_DIRECTORY} DESTINATION ${PROJECT_BINARY_DIR}/resources)
    endif(BUILD_FOR_WEB_DEPLOYMENT)
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

# TODO: I think we can get rid of this next part, but need to modify in the code a little bit probably
if(INCLUDE_ANALYSIS_TEST_SUITE)
    if(WIN32)
        #        file(COPY ${PROJECT_SOURCE_DIR}/target/testing/analysis_tests DESTINATION ${PROJECT_BINARY_DIR}/analysis_tests )
    else(WIN32)
        execute_process(
            COMMAND
                ${CMAKE_COMMAND} -E create_symlink
                ${PROJECT_SOURCE_DIR}/target/testing/analysis_tests
                ${PROJECT_BINARY_DIR}/analysis_tests
        )
    endif(WIN32)
endif()

configure_file(
    ${PROJECT_SOURCE_DIR}//external_libs/SandiaDecay/sandia.decay.min.xml
    ${PROJECT_BINARY_DIR}/data/sandia.decay.xml
    COPYONLY
)

configure_file(
    ${PROJECT_SOURCE_DIR}/data/xray_widths.xml
    ${PROJECT_BINARY_DIR}/data/xray_widths.xml
    COPYONLY
)

# Only create an install target if this is a top-level project (e.g., not Electron, or iOS, etc)
if( NOT hasParent )
  INSTALL(TARGETS InterSpecExe DESTINATION bin)
  
  if( INTERSPEC_LIB_TYPE STREQUAL "SHARED" )
    INSTALL(TARGETS InterSpecLib DESTINATION lib)
  endif()
  
  INSTALL(DIRECTORY data example_spectra InterSpec_resources DESTINATION share/interspec)
  INSTALL(DIRECTORY ${WT_RESOURCES_DIRECTORY}/ DESTINATION share/interspec/resources)
endif( NOT hasParent )
