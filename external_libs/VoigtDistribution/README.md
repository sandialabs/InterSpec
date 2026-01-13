# VoigtDistribution

A C++ library for Voigt profile distributions with exponential low-energy tail, for gamma spectroscopy.  Computation is templated on the computation type so you can use it with either `double`, or `ceres::Jet<>` (e.g. for use with Ceres Solver, http://ceres-solver.org) types.

## Components

- **Faddeeva functions**: High-precision implementation of the Faddeeva function w(z) and related complex error functions - the implementation is a strate-forwardly templated version of the MIT licensed code from Steven G. Johnson, avaiable at http://ab-initio.mit.edu/faddeeva/; the code was last updated 12 May 2015.
- **Voigt profiles**: True Voigt distributions (Gaussian-Lorentzian convolution) using Faddeeva functions
- **Exponential low-energy tail**: Exponentially-modified Gaussian distributions for modeling incomplete charge collection
- **Combined distributions**: Voigt profiles with exponential low-energy tails

## Features

- Template-based implementation supporting `double` and `ceres::Jet<>` for automatic differentiation
- Unit-area normalization
- High numerical accuracy
- Self-contained with minimal dependencies (only standard C++ libraries and optionally Ceres/Eigen), except to run unit tests, which use the Boost unit_test framework.

## Usage

```cpp
#include "VoigtDistribution/voigt_exp_tail.hpp"

// Create a Voigt distribution with exponential tail
double pdf_val = voigt_exp_indefinite(x, mean, sigma, gamma, tail_ratio, tail_slope);
```

## Building

This library is header-only and requires C++17. For Ceres Jet support, link against Ceres and Eigen.