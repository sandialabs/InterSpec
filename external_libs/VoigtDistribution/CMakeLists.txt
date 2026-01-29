cmake_minimum_required(VERSION 3.10)

project(VoigtDistribution)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to build unit tests
option(BUILD_VOIGT_TESTS "Build unit tests for VoigtDistribution" OFF)

# Option to build example
option(BUILD_VOIGT_EXAMPLE "Build example program for VoigtDistribution" OFF)

# Create interface library (header-only)
add_library(VoigtDistribution INTERFACE)
target_include_directories(VoigtDistribution INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

if(BUILD_VOIGT_TESTS)
    # Add unit tests subdirectory
    add_subdirectory(unit_tests)
endif()

if(BUILD_VOIGT_EXAMPLE)
    # Add example executable
    add_executable(voigt_example example/example.cpp)
    target_link_libraries(voigt_example VoigtDistribution)
    target_include_directories(voigt_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()

# Installation (optional)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
    PATTERN "unit_tests" EXCLUDE
    PATTERN "CMakeLists.txt" EXCLUDE
)

# Install headers with proper directory structure
install(FILES voigt_exp_tail.hpp pseudo_voigt_exp_tail.hpp Faddeeva.hpp Faddeeva_impl.hpp
    DESTINATION include/VoigtDistribution
)

install(TARGETS VoigtDistribution
    EXPORT VoigtDistributionTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)