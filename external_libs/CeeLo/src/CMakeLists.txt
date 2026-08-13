add_library(CeeLo
    cross_sections/CrossSectionData.cpp
    cross_sections/element_data.cpp
    cross_sections/photon_epics_data.cpp
    cross_sections/relaxation_epics_data.cpp
    materials/Material.cpp
    materials/BuiltinMaterials.cpp
    geometry/Cylinder.cpp
    geometry/Box.cpp
    geometry/RayTrace.cpp
    geometry/SourceGeometry.cpp
    transport/ComptonScatter.cpp
    transport/PhotonTransport.cpp
    transport/TransportUtils.cpp
    physics/ElectronCsda.cpp
    physics/estar_stopping_data.cpp
    cascade/AnalyticCascade.cpp
    cascade/BetaSpectrum.cpp
    efficiency/EfficiencyCalculator.cpp
    efficiency/SourceDimensionEstimator.cpp
    export/Geant4Export.cpp
    io/VirtualDepthFit.cpp
)

# Storable detector-response payload (optional; see CEELO_BUILD_RESPONSE in the
# top-level CMakeLists, which also locates rapidxml in CEELO_RAPIDXML_DIR).
if(CEELO_BUILD_RESPONSE)
    target_sources(CeeLo PRIVATE
        io/DetectorResponse.cpp
        io/EfficiencyTransfer.cpp
        io/Pchip.cpp
        io/ResponseGenerator.cpp
        io/ResponseKernel.cpp
    )
    target_include_directories(CeeLo PRIVATE ${CEELO_RAPIDXML_DIR})
endif()

target_include_directories(CeeLo
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# The library itself includes no Boost headers -- Boost is a test-only
# dependency (Boost.Test), found in tests/CMakeLists.txt.
target_link_libraries(CeeLo
    PUBLIC
        Eigen3::Eigen
)

# Strict warnings
target_compile_options(CeeLo PRIVATE
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# Optional SandiaDecay -> DecayCascade adapter (CMake-gated; keeps the core
# library free of any SandiaDecay dependency). Built as a separate static
# library that links both CeeLo and SandiaDecay.
if(CEELO_WITH_SANDIADECAY)
    add_library(CeeLoCascade STATIC cascade/SandiaDecayCascade.cpp)
    target_link_libraries(CeeLoCascade PUBLIC CeeLo SandiaDecay)
    target_compile_definitions(CeeLoCascade PUBLIC CEELO_WITH_SANDIADECAY)
    # Not -Wpedantic here: SandiaDecay.h is a third-party header.
    target_compile_options(CeeLoCascade PRIVATE
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall -Wextra>
        $<$<CXX_COMPILER_ID:MSVC>:/W4>
    )
endif()
