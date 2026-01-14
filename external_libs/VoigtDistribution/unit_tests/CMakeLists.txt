cmake_minimum_required(VERSION 3.10)

project(VoigtDistributionUnitTests)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find dependencies needed for unit tests
find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# Enable testing
enable_testing()

# Create unit test executable
add_executable(unit_tests
    voigt_distribution_tests.cpp
    test_faddeeva_templated.cpp
    test_voigt_exact.cpp
    test_pseudo_voigt_exp_tail.cpp
    exact_test_data.cpp
    3rd_party/Faddeeva_MIT/Faddeeva.cc
)

target_include_directories(unit_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party
    ${EIGEN3_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(unit_tests
    Boost::unit_test_framework
    Eigen3::Eigen
    Ceres::ceres
)

# Add test
add_test(NAME VoigtDistributionUnitTests COMMAND unit_tests)